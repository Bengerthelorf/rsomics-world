use std::path::PathBuf;
use std::process::Command;
fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-makewindows"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}
#[test]
fn tiles_100bp_genome_into_25bp_windows() {
    let out = Command::new(bin())
        .args(["-g"])
        .arg(fixture("tiny.genome"))
        .args(["-w", "25"])
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    let count = s.trim().lines().count();
    assert_eq!(count, 4, "100bp / 25bp = 4 windows: {s}");
}
