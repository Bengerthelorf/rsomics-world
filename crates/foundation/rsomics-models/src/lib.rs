#![allow(clippy::missing_errors_doc, clippy::missing_panics_doc)]

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
#[non_exhaustive]
pub enum Format {
    Onnx,
    Safetensors,
    Pt,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Model {
    pub name: String,
    pub version: String,
    pub url: String,
    pub sha256: String,
    pub format: Format,
    /// Cache path stem (without extension); the loader appends the format-
    /// appropriate suffix.
    pub stem: String,
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ModelError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("checksum mismatch for {name}: expected {expected}, got {actual}")]
    ChecksumMismatch {
        name: String,
        expected: String,
        actual: String,
    },
    #[error("model not in cache: {name} (expected at {path})")]
    NotInCache { name: String, path: String },
    #[error("malformed sha256 (expected 64 hex chars, got {0})")]
    BadSha256(String),
}

pub type Result<T> = std::result::Result<T, ModelError>;

#[derive(Debug, Clone)]
pub struct Cache {
    pub root: PathBuf,
}

impl Cache {
    /// Platform cache root, e.g. `~/.cache/rsomics/models` (Linux) or
    /// `~/Library/Caches/rsomics/models` (macOS).
    #[must_use]
    pub fn default_root() -> PathBuf {
        if let Some(home) = std::env::var_os("HOME") {
            let mut p = PathBuf::from(home);
            #[cfg(target_os = "macos")]
            p.push("Library/Caches/rsomics/models");
            #[cfg(not(target_os = "macos"))]
            p.push(".cache/rsomics/models");
            p
        } else {
            PathBuf::from("/tmp/rsomics-models")
        }
    }

    #[must_use]
    pub fn new(root: PathBuf) -> Self {
        Self { root }
    }

    pub fn ensure(&self) -> Result<()> {
        std::fs::create_dir_all(&self.root)?;
        Ok(())
    }

    /// Absolute path to `model` in this cache (present or not).
    #[must_use]
    pub fn path_for(&self, model: &Model) -> PathBuf {
        let ext = match model.format {
            Format::Onnx => "onnx",
            Format::Safetensors => "safetensors",
            Format::Pt => "pt",
        };
        self.root.join(format!("{}.{ext}", model.stem))
    }

    /// `Some(path)` if cached and sha256 matches; `None` if absent;
    /// `Err(ChecksumMismatch)` if present but corrupted.
    pub fn lookup(&self, model: &Model) -> Result<Option<PathBuf>> {
        let p = self.path_for(model);
        if !p.exists() {
            return Ok(None);
        }
        validate_sha256(model)?;
        let actual = sha256_of(&p)?;
        if actual.eq_ignore_ascii_case(&model.sha256) {
            Ok(Some(p))
        } else {
            Err(ModelError::ChecksumMismatch {
                name: model.name.clone(),
                expected: model.sha256.clone(),
                actual,
            })
        }
    }

    /// Like `lookup` but errors if the model is absent.
    pub fn require(&self, model: &Model) -> Result<PathBuf> {
        match self.lookup(model)? {
            Some(p) => Ok(p),
            None => Err(ModelError::NotInCache {
                name: model.name.clone(),
                path: self.path_for(model).display().to_string(),
            }),
        }
    }
}

fn validate_sha256(model: &Model) -> Result<()> {
    if model.sha256.len() != 64 || !model.sha256.chars().all(|c| c.is_ascii_hexdigit()) {
        return Err(ModelError::BadSha256(model.sha256.clone()));
    }
    Ok(())
}

fn sha256_of(path: &Path) -> Result<String> {
    let mut hasher = Sha256::new();
    let mut buf = vec![0_u8; 64 * 1024];
    let mut reader = BufReader::new(File::open(path)?);
    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(hex_encode(&hasher.finalize()))
}

fn hex_encode(bytes: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut out = String::with_capacity(bytes.len() * 2);
    for &b in bytes {
        out.push(HEX[(b >> 4) as usize] as char);
        out.push(HEX[(b & 0xf) as usize] as char);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn make_model(stem: &str, sha: &str) -> Model {
        Model {
            name: stem.to_string(),
            version: "0.1.0".to_string(),
            url: "https://example.org/model.onnx".to_string(),
            sha256: sha.to_string(),
            format: Format::Onnx,
            stem: stem.to_string(),
        }
    }

    fn sha256_of_bytes(bytes: &[u8]) -> String {
        let mut h = Sha256::new();
        h.update(bytes);
        hex_encode(&h.finalize())
    }

    #[test]
    fn lookup_returns_none_when_absent() {
        let dir = tempfile::tempdir().unwrap();
        let cache = Cache::new(dir.path().to_path_buf());
        cache.ensure().unwrap();
        let m = make_model("absent", &"0".repeat(64));
        assert!(cache.lookup(&m).unwrap().is_none());
    }

    #[test]
    fn lookup_validates_existing_file() {
        let dir = tempfile::tempdir().unwrap();
        let cache = Cache::new(dir.path().to_path_buf());
        cache.ensure().unwrap();
        let body = b"hello world model";
        let m = make_model("toy", &sha256_of_bytes(body));
        let p = cache.path_for(&m);
        std::fs::write(&p, body).unwrap();
        let got = cache.lookup(&m).unwrap();
        assert_eq!(got, Some(p));
    }

    #[test]
    fn checksum_mismatch_surfaces_loud() {
        let dir = tempfile::tempdir().unwrap();
        let cache = Cache::new(dir.path().to_path_buf());
        cache.ensure().unwrap();
        let m = make_model("corrupted", &"f".repeat(64));
        let p = cache.path_for(&m);
        std::fs::write(&p, b"not the body the checksum belongs to").unwrap();
        let err = cache.lookup(&m).unwrap_err();
        assert!(matches!(err, ModelError::ChecksumMismatch { .. }));
    }

    #[test]
    fn require_errors_when_missing() {
        let dir = tempfile::tempdir().unwrap();
        let cache = Cache::new(dir.path().to_path_buf());
        cache.ensure().unwrap();
        let m = make_model("missing", &"0".repeat(64));
        assert!(matches!(
            cache.require(&m),
            Err(ModelError::NotInCache { .. })
        ));
    }

    #[test]
    fn bad_sha256_rejected_before_io() {
        let dir = tempfile::tempdir().unwrap();
        let cache = Cache::new(dir.path().to_path_buf());
        cache.ensure().unwrap();
        let m = make_model("bad", "not-64-hex-chars");
        let p = cache.path_for(&m);
        let mut f = File::create(&p).unwrap();
        f.write_all(b"anything").unwrap();
        assert!(matches!(cache.lookup(&m), Err(ModelError::BadSha256(_))));
    }

    #[test]
    fn extension_matches_format() {
        let dir = tempfile::tempdir().unwrap();
        let cache = Cache::new(dir.path().to_path_buf());
        let mut m = make_model("x", &"0".repeat(64));
        m.format = Format::Safetensors;
        assert_eq!(
            cache.path_for(&m).extension().and_then(|s| s.to_str()),
            Some("safetensors")
        );
        m.format = Format::Pt;
        assert_eq!(
            cache.path_for(&m).extension().and_then(|s| s.to_str()),
            Some("pt")
        );
    }
}
