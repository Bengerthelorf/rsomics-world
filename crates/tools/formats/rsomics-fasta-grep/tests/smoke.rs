use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-grep"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn grep_chr1() {
    let out = Command::new(bin())
        .args(["-p", "^chr1"])
        .arg(fixture("three.fa"))
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    let count = s.lines().filter(|l| l.starts_with('>')).count();
    assert_eq!(count, 2, "chr1 matches 2 seqs: {s}");
}

#[test]
fn grep_invert() {
    let out = Command::new(bin())
        .args(["-p", "^chr1", "--invert-match"])
        .arg(fixture("three.fa"))
        .output()
        .expect("spawn");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let count = s.lines().filter(|l| l.starts_with('>')).count();
    assert_eq!(count, 1, "invert: 1 non-chr1 seq: {s}");
}
