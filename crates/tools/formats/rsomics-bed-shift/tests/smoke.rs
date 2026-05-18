use std::path::PathBuf;
use std::process::Command;
fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-shift"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}
#[test]
fn shift_right_50() {
    let out = Command::new(bin())
        .args(["-i"])
        .arg(fixture("one.bed"))
        .args(["-g"])
        .arg(fixture("genome.txt"))
        .args(["-s", "50"])
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    assert_eq!(s.trim(), "chr1\t150\t250", "100+50=150, 200+50=250: {s}");
}
