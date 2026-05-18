use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf { PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-filter")) }
fn fixture(name: &str) -> PathBuf { PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden").join(name) }

#[test]
fn filter_by_min_length() {
    let out = Command::new(bin())
        .args(["-m", "20"])
        .arg(fixture("mixed.fa"))
        .output().expect("spawn");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let count = s.lines().filter(|l| l.starts_with('>')).count();
    assert_eq!(count, 5, "5 long seqs >= 20bp: got {count}");
}

#[test]
fn filter_keeps_all_when_min_zero() {
    let out = Command::new(bin())
        .args(["-m", "0"])
        .arg(fixture("mixed.fa"))
        .output().expect("spawn");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let count = s.lines().filter(|l| l.starts_with('>')).count();
    assert_eq!(count, 10, "min=0 keeps all 10: got {count}");
}
