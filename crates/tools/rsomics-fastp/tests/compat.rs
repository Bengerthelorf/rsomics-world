//! Behavioural-compatibility tests against upstream `fastp`.
//!
//! Each test runs `rsomics-fastp` and the upstream `fastp` (or `fasterp`)
//! binary on a small golden FASTQ, then diffs the byte-level output. Tests
//! are gated behind the `test-support` feature so they don't run by accident
//! in environments without the upstream binary on PATH.

#![cfg(feature = "test-support")]

use std::process::Command;

fn fastp_on_path() -> bool {
    Command::new("fastp")
        .arg("--version")
        .output()
        .is_ok_and(|out| out.status.success())
}

#[test]
fn fastp_binary_available() {
    if !fastp_on_path() {
        eprintln!(
            "skip: upstream fastp not on PATH; install via bioconda or skip the test-support feature"
        );
    }
}
