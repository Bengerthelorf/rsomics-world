use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-to-fastq"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn converts_with_quality() {
    let out = Command::new(bin())
        .arg(fixture("one.fa"))
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    assert!(s.starts_with("@seq1\n"), "FASTQ header: {s}");
    assert!(s.contains("ACGTACGT"), "sequence preserved");
    assert!(s.contains("+\n"), "quality separator");
    let qual_line = s.lines().nth(3).unwrap();
    assert_eq!(qual_line.len(), 8, "quality same length as seq");
}
