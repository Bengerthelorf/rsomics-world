use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-len"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn correct_lengths() {
    let out = Command::new(bin())
        .arg(fixture("two.fq"))
        .output()
        .expect("spawn");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let lengths: Vec<&str> = s.trim().lines().collect();
    assert_eq!(lengths, vec!["4", "10"]);
}
