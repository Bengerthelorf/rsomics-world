use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-closest"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn finds_closest() {
    let out = Command::new(bin())
        .args(["-a"])
        .arg(fixture("a.bed"))
        .args(["-b"])
        .arg(fixture("b.bed"))
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    assert!(
        s.contains("300"),
        "closest to 100-200 should be 300-400: {s}"
    );
}
