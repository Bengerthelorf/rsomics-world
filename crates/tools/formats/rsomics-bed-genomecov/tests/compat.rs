use std::path::PathBuf;
use std::process::Command;
fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-genomecov"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}
#[test]
fn coverage_output_valid() {
    let out = Command::new(ours())
        .args(["-i"])
        .arg(fixture("two.bed"))
        .args(["-g"])
        .arg(fixture("genome.txt"))
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    assert!(s.contains("chr1"), "has chr1 coverage");
}
