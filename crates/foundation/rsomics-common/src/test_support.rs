//! Test-only scaffolding shared across `rsomics-*` integration tests and
//! benches. Behind the `test-support` Cargo feature so the production binary
//! never carries the symbols.
//!
//! The fixture-path resolver is a macro rather than a function because
//! `env!("CARGO_MANIFEST_DIR")` resolves at the *call site* — we need each
//! consuming crate to read its OWN manifest dir, not rsomics-common's.

use std::process::{Command, Output, Stdio};

/// Build a `PathBuf` pointing to `tests/golden/<rel>` under the **caller's**
/// crate manifest directory.
///
/// ```no_run
/// # use rsomics_common::fixture_path;
/// let p = fixture_path!("tiny.bam");
/// assert!(p.ends_with("tests/golden/tiny.bam"));
/// ```
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

/// True if `<name>` is on `PATH` and `<name> --version` exits zero. Used by
/// compat tests to skip gracefully when the upstream binary isn't installed.
#[must_use]
pub fn tool_on_path(name: &str) -> bool {
    Command::new(name)
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

/// Run an external tool with the given args and return its captured `Output`.
/// Panics if the spawn fails entirely — callers should have gated this with
/// [`tool_on_path`] first.
///
/// # Panics
///
/// Panics if `Command::output` returns `Err` (the binary couldn't even be
/// spawned — distinct from "spawned and exited non-zero").
#[must_use]
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
        assert!(tool_on_path("cargo"), "cargo must be on PATH for cargo test");
    }

    #[test]
    fn tool_on_path_returns_false_for_unlikely_binary() {
        assert!(!tool_on_path("rsomics-deliberately-not-installed-xyzzy"));
    }
}
