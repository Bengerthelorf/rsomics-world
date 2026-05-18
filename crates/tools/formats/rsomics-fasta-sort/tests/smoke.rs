use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-sort"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn sort_by_name() {
    let out = Command::new(bin())
        .arg(fixture("unsorted.fa"))
        .output()
        .expect("spawn");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let names: Vec<&str> = s.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(names, vec![">aaa", ">mmm", ">zzz"]);
}

#[test]
fn sort_by_length_desc() {
    let out = Command::new(bin())
        .args(["-L"])
        .arg(fixture("unsorted.fa"))
        .output()
        .expect("spawn");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let first_name = s.lines().next().unwrap();
    assert_eq!(
        first_name, ">aaa",
        "longest (8bp) should be first: got {first_name}"
    );
}
