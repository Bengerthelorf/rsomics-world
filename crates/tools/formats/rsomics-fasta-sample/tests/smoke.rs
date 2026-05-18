use std::path::PathBuf;
use std::process::Command;
fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-sample"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}
#[test]
fn samples_subset() {
    let out = Command::new(bin())
        .args(["-p", "0.5", "--seed", "42"])
        .arg(fixture("ten.fa"))
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
        .filter(|l| l.starts_with('>'))
        .count();
    assert!(count > 0 && count < 10, "subset of 10: got {count}");
}
