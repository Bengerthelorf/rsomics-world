#[cfg(feature = "tier2")]
pub mod tier2;

use std::process::{Command, Output, Stdio};

/// `PathBuf` pointing to `tests/golden/<rel>` under the caller's crate.
/// Macro so `env!("CARGO_MANIFEST_DIR")` resolves at the call site.
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
        p.pop(); // crates/{foundation,tools}/<crate> → crates/{foundation,tools}
        p.pop(); // → crates
        p.pop(); // → workspace root
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
