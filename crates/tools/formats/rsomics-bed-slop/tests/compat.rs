use std::path::PathBuf;
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-slop"))
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
fn slop_matches_bedtools() {
    if !bedtools_available() {
        eprintln!("bedtools not on PATH — skipping compat test");
        return;
    }
    let bed = fixture("two.bed");
    let genome = fixture("genome.txt");

    let ours_out = Command::new(ours())
        .args(["-i"])
        .arg(&bed)
        .args(["-g"])
        .arg(&genome)
        .args(["-b", "10"])
        .output()
        .expect("spawn ours");
    assert!(ours_out.status.success());

    let bt_out = Command::new("bedtools")
        .args(["slop", "-i"])
        .arg(&bed)
        .args(["-g"])
        .arg(&genome)
        .args(["-b", "10"])
        .output()
        .expect("spawn bedtools");
    assert!(bt_out.status.success());

    let ours_str = String::from_utf8(ours_out.stdout).unwrap();
    let bt_str = String::from_utf8(bt_out.stdout).unwrap();

    assert_eq!(
        ours_str.trim(),
        bt_str.trim(),
        "output must match bedtools:\nours:\n{ours_str}\nbedtools:\n{bt_str}"
    );
}
