//! Tier-2 fixture management: manifest-driven download + sha256 verification.
//! Behind the `tier2` Cargo feature so production binaries don't carry
//! `ureq` / `sha2` / `toml`.

use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use serde::Deserialize;
use sha2::{Digest, Sha256};

use crate::error::{Context, Result, RsomicsError};

#[derive(Debug, Clone, Deserialize)]
pub struct Fixture {
    pub name: String,
    pub tier: u8,
    pub url: String,
    pub sha256: String,
    pub size_mb: u64,
    pub license: String,
    pub source: String,
}

#[derive(Debug, Deserialize)]
struct Manifest {
    #[serde(default, rename = "fixture")]
    fixtures: Vec<Fixture>,
}

/// Resolve the Tier-2 cache directory, creating it if absent.
///
/// Precedence: `RSOMICS_TIER2_CACHE` > `$CARGO_TARGET_DIR/tier2-cache` > `~/.cache/rsomics-fixtures`.
///
/// # Errors
///
/// `RsomicsError::Io` if mkdir fails; `RsomicsError::ConfigError` if `HOME` is unset.
pub fn cache_dir() -> Result<PathBuf> {
    resolve_cache_dir(
        std::env::var("RSOMICS_TIER2_CACHE").ok().as_deref(),
        std::env::var("CARGO_TARGET_DIR").ok().as_deref(),
        std::env::var("HOME").ok().as_deref(),
    )
}

fn resolve_cache_dir(
    override_dir: Option<&str>,
    cargo_target: Option<&str>,
    home: Option<&str>,
) -> Result<PathBuf> {
    let dir = if let Some(p) = override_dir {
        PathBuf::from(p)
    } else if let Some(t) = cargo_target {
        PathBuf::from(t).join("tier2-cache")
    } else {
        let home = home.ok_or_else(|| {
            RsomicsError::ConfigError("HOME unset; cannot resolve cache dir".into())
        })?;
        PathBuf::from(home).join(".cache").join("rsomics-fixtures")
    };
    std::fs::create_dir_all(&dir)
        .rs_with_context(|| format!("creating tier-2 cache dir {}", dir.display()))?;
    Ok(dir)
}

/// # Errors
///
/// `Err` if the file can't be read or doesn't parse.
pub fn load_manifest(path: &Path) -> Result<Vec<Fixture>> {
    let mut buf = String::new();
    File::open(path)
        .rs_with_context(|| format!("opening manifest {}", path.display()))?
        .read_to_string(&mut buf)
        .rs_with_context(|| format!("reading manifest {}", path.display()))?;
    let m: Manifest = toml::from_str(&buf)
        .map_err(|e| RsomicsError::InvalidInput(format!("parsing manifest: {e}")))?;
    Ok(m.fixtures)
}

/// Ensure fixture `name` is present and sha256-verified in the cache, downloading if needed.
/// Verification happens before the rename — a partial download cannot poison the cache.
///
/// # Errors
///
/// `Err` if not found in manifest, download fails, sha256 mismatches, or filesystem ops fail.
pub fn fetch(manifest_path: &Path, name: &str) -> Result<PathBuf> {
    let fixtures = load_manifest(manifest_path)?;
    let fx = fixtures.iter().find(|f| f.name == name).ok_or_else(|| {
        RsomicsError::InvalidInput(format!("fixture {name:?} not found in manifest"))
    })?;
    let dest = cache_dir()?.join(&fx.name);

    if dest.exists()
        && let Ok(actual) = sha256_of(&dest)
        && actual.eq_ignore_ascii_case(&fx.sha256)
    {
        return Ok(dest);
    }

    let tmp = dest.with_extension("part");
    download_to(&fx.url, &tmp)?;
    let got = sha256_of(&tmp)?;
    if !got.eq_ignore_ascii_case(&fx.sha256) {
        let _ = std::fs::remove_file(&tmp);
        return Err(RsomicsError::InvalidInput(format!(
            "sha256 mismatch for {name}: expected {}, got {got}",
            fx.sha256
        )));
    }
    std::fs::rename(&tmp, &dest)
        .rs_with_context(|| format!("moving fixture into cache: {}", dest.display()))?;
    Ok(dest)
}

fn sha256_of(path: &Path) -> Result<String> {
    let mut f = File::open(path).rs_with_context(|| format!("opening {}", path.display()))?;
    let mut hasher = Sha256::new();
    let mut buf = vec![0u8; 64 * 1024];
    loop {
        let n = f
            .read(&mut buf)
            .rs_with_context(|| format!("reading {}", path.display()))?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

fn download_to(url: &str, dest: &Path) -> Result<()> {
    let mut response = ureq::get(url)
        .call()
        .map_err(|e| RsomicsError::UpstreamError(format!("HTTP GET {url}: {e}")))?;
    let mut reader = response.body_mut().as_reader();
    let mut out = File::create(dest)
        .rs_with_context(|| format!("creating partial download {}", dest.display()))?;
    let mut buf = vec![0u8; 128 * 1024];
    loop {
        let n = reader
            .read(&mut buf)
            .map_err(|e| RsomicsError::UpstreamError(format!("HTTP read body for {url}: {e}")))?;
        if n == 0 {
            break;
        }
        out.write_all(&buf[..n])
            .rs_with_context(|| format!("writing {}", dest.display()))?;
    }
    out.flush()
        .rs_with_context(|| format!("flushing {}", dest.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sha256_matches_known_vector() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("hello.txt");
        std::fs::write(&path, b"hello\n").expect("write");
        // sha256("hello\n") = 5891b5b522d5df086d0ff0b110fbd9d21bb4fc7163af34d08286a2e846f6be03
        let got = sha256_of(&path).expect("hash");
        assert_eq!(
            got,
            "5891b5b522d5df086d0ff0b110fbd9d21bb4fc7163af34d08286a2e846f6be03"
        );
    }

    #[test]
    fn manifest_parses_round_trip() {
        let dir = tempfile::tempdir().expect("tempdir");
        let manifest = dir.path().join("manifest.toml");
        std::fs::write(
            &manifest,
            r#"
[[fixture]]
name = "hg002-chr22-fastq"
tier = 2
url = "https://example.com/x.fastq.gz"
sha256 = "deadbeef"
size_mb = 100
license = "CC0-1.0"
source = "GIAB"
"#,
        )
        .expect("write");
        let fxs = load_manifest(&manifest).expect("load");
        assert_eq!(fxs.len(), 1);
        assert_eq!(fxs[0].name, "hg002-chr22-fastq");
        assert_eq!(fxs[0].tier, 2);
    }

    #[test]
    fn manifest_unknown_fixture_errors() {
        let dir = tempfile::tempdir().expect("tempdir");
        let manifest = dir.path().join("manifest.toml");
        std::fs::write(&manifest, "").expect("write");
        let r = fetch(&manifest, "nope");
        assert!(r.is_err());
    }

    #[test]
    fn resolve_uses_explicit_override_first() {
        let dir = tempfile::tempdir().expect("tempdir");
        let override_path = dir.path().to_string_lossy().into_owned();
        let got = resolve_cache_dir(
            Some(&override_path),
            Some("/no/such/target"),
            Some("/no/home"),
        )
        .expect("resolve");
        assert_eq!(got, dir.path());
        assert!(got.is_dir());
    }

    #[test]
    fn resolve_falls_back_to_cargo_target_dir() {
        let dir = tempfile::tempdir().expect("tempdir");
        let target = dir.path().to_string_lossy().into_owned();
        let got = resolve_cache_dir(None, Some(&target), Some("/no/home")).expect("resolve");
        assert_eq!(got, dir.path().join("tier2-cache"));
        assert!(got.is_dir());
    }

    #[test]
    fn resolve_falls_back_to_home_cache() {
        let dir = tempfile::tempdir().expect("tempdir");
        let home = dir.path().to_string_lossy().into_owned();
        let got = resolve_cache_dir(None, None, Some(&home)).expect("resolve");
        assert_eq!(got, dir.path().join(".cache").join("rsomics-fixtures"));
        assert!(got.is_dir());
    }

    #[test]
    fn resolve_errors_when_all_unset() {
        let r = resolve_cache_dir(None, None, None);
        assert!(matches!(r, Err(RsomicsError::ConfigError(_))));
    }
}
