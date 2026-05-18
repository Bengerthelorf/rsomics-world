use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf { PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-len")) }
fn fixture(name: &str) -> PathBuf { PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden").join(name) }

#[test]
fn outputs_correct_lengths() {
    let out = Command::new(bin()).arg(fixture("three.fa")).output().expect("spawn");
    assert!(out.status.success());
    let lengths: Vec<&str> = String::from_utf8(out.stdout).unwrap().trim().split('\n').collect::<Vec<_>>().into_iter().collect();
    assert_eq!(lengths, vec!["4", "8", "2"]);
}

#[test]
fn tab_mode_includes_names() {
    let out = Command::new(bin()).args(["--tab"]).arg(fixture("three.fa")).output().expect("spawn");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    assert!(s.contains("a\t4"), "tab mode: {s}");
    assert!(s.contains("b\t8"), "tab mode: {s}");
}
