use std::path::PathBuf;
use std::process::Command;
fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-split"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn splits_into_multiple_files() {
    let tmp = tempfile::tempdir().unwrap();
    let prefix = format!("{}/part_", tmp.path().display());
    let out = Command::new(bin())
        .args(["--seqs-per-file", "2", "--prefix", &prefix])
        .arg(fixture("five.fa"))
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let files: Vec<_> = std::fs::read_dir(tmp.path()).unwrap().collect();
    assert_eq!(files.len(), 3, "5 seqs / 2 per file = 3 files");
}
