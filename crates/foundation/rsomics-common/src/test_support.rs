
#[cfg(feature = "tier2")]
pub mod tier2;

use std::process::{Command, Output, Stdio};

/// Build a `PathBuf` pointing to `tests/golden/<rel>` under the **caller's**
/// crate manifest directory. Must be a macro so `env!("CARGO_MANIFEST_DIR")`
/// resolves at the call site.
#[macro_export]
macro_rules! fixture_path {
    ($rel:expr) => {{
        let mut p = ::std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        p.push("tests");
        p.push("golden");
        p.push($rel);
        p
    }};
}

#[macro_export]
macro_rules! tier2_manifest_path {
    () => {{
        let mut p = ::std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        // CARGO_MANIFEST_DIR = <root>/crates/{foundation,tools}/<crate>
        p.pop(); // <root>/crates/{foundation,tools}
        p.pop(); // <root>/crates
        p.pop(); // <root>
        p.push("tests");
        p.push("fixtures-manifest.toml");
        p
    }};
}

#[must_use]
pub fn tool_on_path(name: &str) -> bool {
    Command::new(name)
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

#[must_use]
#[allow(clippy::missing_panics_doc)]
pub fn run_tool(name: &str, args: &[&str]) -> Output {
    Command::new(name)
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("spawn {name}: {e}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fixture_path_macro_resolves_under_manifest_dir() {
        let p: std::path::PathBuf = fixture_path!("nonexistent.txt");
        let s = p.to_string_lossy();
        assert!(s.ends_with("tests/golden/nonexistent.txt"), "got {s}");
        assert!(s.contains("rsomics-common"), "got {s}");
    }

    #[test]
    fn tool_on_path_finds_cargo() {
        // `cargo` is universally present during `cargo test` and supports
        // `--version`. POSIX `sh` does not — on Debian/Ubuntu `sh` is
        // `dash`, which exits non-zero on `--version`. The predicate is
        // intentionally aimed at tools that follow the `--version`
        // convention (samtools, fastp, bcftools, …); using a shell that
        // doesn't would falsely report it absent.
        assert!(
            tool_on_path("cargo"),
            "cargo must be on PATH for cargo test"
        );
    }

    #[test]
    fn tool_on_path_returns_false_for_unlikely_binary() {
        assert!(!tool_on_path("rsomics-deliberately-not-installed-xyzzy"));
    }
}
