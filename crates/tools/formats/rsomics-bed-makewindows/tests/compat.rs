use std::path::PathBuf;
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-makewindows"))
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
fn windows_match_bedtools() {
    if !bedtools_available() {
        eprintln!("bedtools not on PATH — skipping");
        return;
    }
    let genome = fixture("tiny.genome");
    let ours_out = Command::new(ours())
        .args(["-g"])
        .arg(&genome)
        .args(["-w", "25"])
        .output()
        .expect("spawn");
    assert!(ours_out.status.success());
    let bt_out = Command::new("bedtools")
        .args(["makewindows", "-g"])
        .arg(&genome)
        .args(["-w", "25"])
        .output()
        .expect("spawn");
    assert!(bt_out.status.success());
    assert_eq!(
        String::from_utf8(ours_out.stdout).unwrap().trim(),
        String::from_utf8(bt_out.stdout).unwrap().trim(),
        "must match bedtools makewindows"
    );
}
