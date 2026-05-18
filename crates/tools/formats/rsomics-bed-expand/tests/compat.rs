use std::path::PathBuf;
use std::process::Command;
fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-expand"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}
#[test]
fn passthrough_correct() {
    let ours_out = Command::new(ours())
        .arg(fixture("one.bed"))
        .output()
        .expect("spawn");
    assert!(ours_out.status.success());
    assert!(
        String::from_utf8(ours_out.stdout)
            .unwrap()
            .contains("chr1\t10\t13")
    );
}
