use std::path::PathBuf;
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-closest"))
}
fn bedtools_available() -> bool {
    Command::new("bedtools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn closest_produces_output() {
    if !bedtools_available() {
        eprintln!("bedtools not on PATH — skipping");
        return;
    }
    let a = fixture("a.bed");
    let b = fixture("b.bed");
    let ours_out = Command::new(ours())
        .args(["-a"])
        .arg(&a)
        .args(["-b"])
        .arg(&b)
        .output()
        .expect("spawn");
    assert!(
        ours_out.status.success(),
        "{}",
        String::from_utf8_lossy(&ours_out.stderr)
    );
    let s = String::from_utf8(ours_out.stdout).unwrap();
    assert!(!s.trim().is_empty(), "should produce output");
    assert!(s.contains("chr1"), "should reference chr1");
}
