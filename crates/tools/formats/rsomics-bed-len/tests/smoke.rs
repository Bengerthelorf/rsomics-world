use std::path::PathBuf;
use std::process::Command;
fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-len"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn correct_lengths() {
    let out = Command::new(bin())
        .arg(fixture("three.bed"))
        .output()
        .expect("spawn");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let lens: Vec<&str> = s.trim().lines().collect();
    assert_eq!(
        lens,
        vec!["10", "20", "100"],
        "10-20=10, 30-50=20, 100-200=100"
    );
}
