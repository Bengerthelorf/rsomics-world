use std::path::PathBuf;
use std::process::Command;
fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-head"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}
#[test]
fn head_2() {
    let out = Command::new(bin())
        .args(["-n", "2"])
        .arg(fixture("five.fq"))
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
    assert_eq!(count, 2);
}
