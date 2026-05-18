use std::path::PathBuf;
use std::process::Command;
fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-random"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}
#[test]
fn generates_100_intervals() {
    let out = Command::new(bin())
        .args(["-g"])
        .arg(fixture("genome.txt"))
        .args(["-n", "100", "-l", "500", "--seed", "42"])
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
    assert_eq!(count, 100, "expected 100 intervals: got {count}");
}
