use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-to-fasta"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn converts_to_fasta() {
    let out = Command::new(bin())
        .arg(fixture("one.fq"))
        .output()
        .expect("spawn");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    assert!(
        s.starts_with(">r0\n"),
        "should start with FASTA header: {s}"
    );
    assert!(s.contains("ACGTACGT"), "should contain sequence");
    assert!(!s.contains('+'), "should not contain quality separator");
}
