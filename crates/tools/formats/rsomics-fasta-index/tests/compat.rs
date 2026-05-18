use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-index"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn index_creates_fai() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("test.fa");
    std::fs::copy(fixture("two.fa"), &fa).unwrap();
    let out = Command::new(ours())
        .args(["index"])
        .arg(&fa)
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let fai = fa.with_extension("fa.fai");
    assert!(fai.exists(), ".fai must be created");
}
