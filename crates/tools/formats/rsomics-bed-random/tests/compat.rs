use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-random"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn generates_correct_count() {
    let out = Command::new(ours())
        .args(["-g"])
        .arg(fixture("genome.txt"))
        .args(["-n", "50", "-l", "100"])
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let count = String::from_utf8(out.stdout)
        .unwrap()
        .trim()
        .lines()
        .count();
    assert_eq!(count, 50, "should generate exactly 50 intervals");
}
