use std::path::PathBuf;
use std::process::Command;
fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-stats"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}
#[test]
fn total_bases_matches_awk() {
    let input = fixture("three.bed");
    let ours_out = Command::new(ours()).arg(&input).output().expect("spawn");
    assert!(ours_out.status.success());
    let s = String::from_utf8(ours_out.stdout).unwrap();
    assert!(s.contains("Total bases:\t350"), "awk: 100+200+50=350: {s}");
}
