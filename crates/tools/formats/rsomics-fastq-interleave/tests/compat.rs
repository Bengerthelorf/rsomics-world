use std::path::PathBuf;
use std::process::Command;
fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-interleave"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}
#[test]
fn interleaves() {
    let out = Command::new(ours())
        .args(["-i"])
        .arg(fixture("r1.fq"))
        .args(["-I"])
        .arg(fixture("r2.fq"))
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let count = String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .filter(|l| l.starts_with('@'))
        .count();
    assert_eq!(count, 4, "2 pairs interleaved = 4 records");
}
