use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-len"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn lengths_match_awk() {
    let input = fixture("three.bed");
    let ours_out = Command::new(ours()).arg(&input).output().expect("spawn");
    assert!(ours_out.status.success());
    let awk_out = Command::new("awk")
        .args(["-F\t", "{print $3-$2}"])
        .arg(&input)
        .output()
        .expect("spawn awk");
    assert!(awk_out.status.success());
    assert_eq!(
        String::from_utf8(ours_out.stdout).unwrap().trim(),
        String::from_utf8(awk_out.stdout).unwrap().trim(),
        "lengths must match awk"
    );
}
