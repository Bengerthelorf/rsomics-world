use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf { PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-rename")) }
fn fixture(name: &str) -> PathBuf { PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden").join(name) }

#[test]
fn renames_with_prefix() {
    let out = Command::new(bin())
        .args(["--prefix", "sample_"])
        .arg(fixture("one.fq"))
        .output().expect("spawn");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let s = String::from_utf8(out.stdout).unwrap();
    assert!(s.starts_with("@sample_0\n"), "should rename: {s}");
    assert!(!s.contains("old_name"), "old name should be gone");
}
