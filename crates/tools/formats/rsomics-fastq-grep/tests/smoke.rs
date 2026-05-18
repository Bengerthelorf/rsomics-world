use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-grep"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn grep_matches_pattern() {
    let out = Command::new(bin())
        .args(["-p", "yes$"])
        .arg(fixture("three.fq"))
        .output()
        .expect("spawn");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let count = s.lines().filter(|l| l.starts_with('@')).count();
    assert_eq!(count, 2, "expected 2 matching reads: {s}");
}

#[test]
fn grep_invert_match() {
    let out = Command::new(bin())
        .args(["-p", "yes$", "--invert-match"])
        .arg(fixture("three.fq"))
        .output()
        .expect("spawn");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let count = s.lines().filter(|l| l.starts_with('@')).count();
    assert_eq!(count, 1, "expected 1 non-matching read: {s}");
}
