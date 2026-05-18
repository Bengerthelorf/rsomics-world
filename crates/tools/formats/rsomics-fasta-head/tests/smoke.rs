use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-head"))
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn default_head_outputs_10_records() {
    let out = Command::new(bin())
        .arg(fixture("twenty.fa"))
        .output()
        .expect("spawn");
    assert!(out.status.success(), "failed: {}", String::from_utf8_lossy(&out.stderr));
    let stdout = String::from_utf8(out.stdout).unwrap();
    let record_count = stdout.lines().filter(|l| l.starts_with('>')).count();
    assert_eq!(record_count, 10, "default -n 10 should output 10 records, got {record_count}");
}

#[test]
fn head_n5_outputs_5_records() {
    let out = Command::new(bin())
        .args(["-n", "5"])
        .arg(fixture("twenty.fa"))
        .output()
        .expect("spawn");
    assert!(out.status.success());
    let stdout = String::from_utf8(out.stdout).unwrap();
    let record_count = stdout.lines().filter(|l| l.starts_with('>')).count();
    assert_eq!(record_count, 5);
}

#[test]
fn head_n0_outputs_nothing() {
    let out = Command::new(bin())
        .args(["-n", "0"])
        .arg(fixture("twenty.fa"))
        .output()
        .expect("spawn");
    assert!(out.status.success());
    assert!(out.stdout.is_empty());
}

#[test]
fn head_n100_outputs_all_20() {
    let out = Command::new(bin())
        .args(["-n", "100"])
        .arg(fixture("twenty.fa"))
        .output()
        .expect("spawn");
    assert!(out.status.success());
    let stdout = String::from_utf8(out.stdout).unwrap();
    let record_count = stdout.lines().filter(|l| l.starts_with('>')).count();
    assert_eq!(record_count, 20, "requesting more than available should output all");
}
