use std::path::PathBuf;
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-flank"))
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
fn flank_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("bedtools not on PATH — skipping");
        return;
    }
    let bed = fixture("one.bed");
    let genome = fixture("genome.txt");
    let ours_out = Command::new(ours())
        .args(["-i"])
        .arg(&bed)
        .args(["-g"])
        .arg(&genome)
        .args(["-b", "50"])
        .output()
        .expect("spawn");
    assert!(ours_out.status.success());
    let bt_out = Command::new("bedtools")
        .args(["flank", "-i"])
        .arg(&bed)
        .args(["-g"])
        .arg(&genome)
        .args(["-b", "50"])
        .output()
        .expect("spawn");
    assert!(bt_out.status.success());
    assert_eq!(
        String::from_utf8(ours_out.stdout).unwrap().trim(),
        String::from_utf8(bt_out.stdout).unwrap().trim(),
        "must match bedtools flank"
    );
}
