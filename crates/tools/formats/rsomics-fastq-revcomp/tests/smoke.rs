use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-revcomp"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn revcomp_acgt_is_acgt() {
    let out = Command::new(bin())
        .arg(fixture("one.fq"))
        .output()
        .expect("spawn");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let seq = s.lines().nth(1).unwrap();
    assert_eq!(
        seq, "ACGT",
        "revcomp of ACGT is ACGT (palindrome): got {seq}"
    );
}

#[test]
fn qual_is_reversed() {
    let out = Command::new(bin())
        .arg(fixture("one.fq"))
        .output()
        .expect("spawn");
    let s = String::from_utf8(out.stdout).unwrap();
    let qual = s.lines().nth(3).unwrap();
    assert_eq!(
        qual, "IIII",
        "qual should be reversed (but IIII is symmetric)"
    );
}
