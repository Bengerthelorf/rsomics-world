use std::path::PathBuf;
use std::process::Command;
fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-genomecov"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}
#[test]
fn coverage_output() {
    let out = Command::new(bin())
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
    assert!(!s.is_empty(), "should produce output");
    assert!(s.contains("chr1"), "should contain chr1");
}
