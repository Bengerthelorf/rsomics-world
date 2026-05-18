use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-count"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn count_matches_grep() {
    let input = fixture("three.bed");
    let ours_out = Command::new(ours()).arg(&input).output().expect("spawn");
    assert!(ours_out.status.success());
    let our_count: u64 = String::from_utf8(ours_out.stdout)
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    let grep_out = Command::new("grep")
        .args(["-cv", "^#"])
        .arg(&input)
        .output()
        .expect("spawn grep");
    let grep_count: u64 = String::from_utf8(grep_out.stdout)
        .unwrap()
        .trim()
        .parse()
        .unwrap_or(0);
    assert_eq!(our_count, grep_count, "must match grep -cv '^#'");
}
