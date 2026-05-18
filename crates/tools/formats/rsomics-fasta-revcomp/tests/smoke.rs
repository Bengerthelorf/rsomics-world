use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-revcomp"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn aaaa_becomes_tttt() {
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
    let seq = s.lines().nth(1).unwrap();
    assert_eq!(seq, "TTTT", "revcomp of AAAA is TTTT: got {seq}");
}
