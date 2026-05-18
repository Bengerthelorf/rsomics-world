use std::path::PathBuf;
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-rename"))
}
fn seqkit_available() -> bool {
    Command::new("seqkit")
        .arg("version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn rename_produces_sequential_ids() {
    if !seqkit_available() {
        eprintln!("seqkit not on PATH — skipping");
        return;
    }
    let input = fixture("one.fa");
    let ours_out = Command::new(ours())
        .args(["--prefix", "seq_"])
        .arg(&input)
        .output()
        .expect("spawn");
    assert!(ours_out.status.success());
    let s = String::from_utf8(ours_out.stdout).unwrap();
    assert!(s.starts_with(">seq_0\n"), "sequential IDs: {s}");
}
