//! Byte-level compat tests vs upstream fastp.
//!
//! Each test runs both binaries on the same input with semantically
//! equivalent flag combinations and diffs the resulting FASTQ. fastp's
//! quality / length / polyG filters are explicitly disabled where they
//! aren't the test's subject, so the comparison isolates the trim layer
//! under test.

use std::path::PathBuf;
use std::process::{Command, Stdio};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-trim"))
}

fn fastp_available() -> bool {
    Command::new("fastp")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

fn run_to_path(bin: &std::path::Path, args: &[&str]) {
    let out = Command::new(bin)
        .args(args)
        .output()
        .expect("subprocess spawn");
    assert!(
        out.status.success(),
        "{} {args:?} failed: stderr=\n{}",
        bin.display(),
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn se_adapter_trim_matches_fastp() {
    assert!(
        fastp_available(),
        "compat test requires fastp on PATH (install via `brew install fastp` / `apt install fastp`)"
    );
    let tmp = tempfile::tempdir().unwrap();
    let ours_out = tmp.path().join("ours.fq");
    let theirs_out = tmp.path().join("theirs.fq");
    let fastp_json = tmp.path().join("fastp.json");
    let fastp_html = tmp.path().join("fastp.html");
    let input = fixture("se_adapter.fastq");

    run_to_path(
        &ours(),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            ours_out.to_str().unwrap(),
            "-a",
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "-L",
        ],
    );
    run_to_path(
        std::path::Path::new("fastp"),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            theirs_out.to_str().unwrap(),
            "-a",
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "-Q",
            "-L",
            "-G",
            "-j",
            fastp_json.to_str().unwrap(),
            "-h",
            fastp_html.to_str().unwrap(),
        ],
    );
    let ours_bytes = std::fs::read(&ours_out).unwrap();
    let theirs_bytes = std::fs::read(&theirs_out).unwrap();
    assert_eq!(ours_bytes, theirs_bytes, "FASTQ output differs");
}

#[test]
fn se_polyg_trim_matches_fastp() {
    assert!(
        fastp_available(),
        "compat test requires fastp on PATH (install via `brew install fastp` / `apt install fastp`)"
    );
    let tmp = tempfile::tempdir().unwrap();
    let ours_out = tmp.path().join("ours.fq");
    let theirs_out = tmp.path().join("theirs.fq");
    let fastp_json = tmp.path().join("fastp.json");
    let fastp_html = tmp.path().join("fastp.html");
    let input = fixture("se_polyg.fastq");

    run_to_path(
        &ours(),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            ours_out.to_str().unwrap(),
            "-A",
            "-g",
            "-L",
        ],
    );
    run_to_path(
        std::path::Path::new("fastp"),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            theirs_out.to_str().unwrap(),
            "-A",
            "-g",
            "-Q",
            "-L",
            "-j",
            fastp_json.to_str().unwrap(),
            "-h",
            fastp_html.to_str().unwrap(),
        ],
    );
    let ours_bytes = std::fs::read(&ours_out).unwrap();
    let theirs_bytes = std::fs::read(&theirs_out).unwrap();
    assert_eq!(ours_bytes, theirs_bytes, "polyG output differs");
}

#[test]
fn se_fixed_trim_matches_fastp() {
    assert!(
        fastp_available(),
        "compat test requires fastp on PATH (install via `brew install fastp` / `apt install fastp`)"
    );
    let tmp = tempfile::tempdir().unwrap();
    let ours_out = tmp.path().join("ours.fq");
    let theirs_out = tmp.path().join("theirs.fq");
    let fastp_json = tmp.path().join("fastp.json");
    let fastp_html = tmp.path().join("fastp.html");
    let input = fixture("se_adapter.fastq");

    run_to_path(
        &ours(),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            ours_out.to_str().unwrap(),
            "-A",
            "-f",
            "4",
            "--trim_tail1",
            "2",
            "-L",
        ],
    );
    run_to_path(
        std::path::Path::new("fastp"),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            theirs_out.to_str().unwrap(),
            "-A",
            "-f",
            "4",
            "-t",
            "2",
            "-Q",
            "-L",
            "-G",
            "-j",
            fastp_json.to_str().unwrap(),
            "-h",
            fastp_html.to_str().unwrap(),
        ],
    );
    let ours_bytes = std::fs::read(&ours_out).unwrap();
    let theirs_bytes = std::fs::read(&theirs_out).unwrap();
    assert_eq!(ours_bytes, theirs_bytes, "fixed-trim output differs");
}
