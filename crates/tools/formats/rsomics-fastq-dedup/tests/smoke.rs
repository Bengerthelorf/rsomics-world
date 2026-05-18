use std::path::PathBuf;
use std::process::Command;
fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-dedup"))
}
#[test]
fn help_exits_zero() {
    let out = Command::new(bin())
        .args(["--help"])
        .output()
        .expect("spawn");
    assert!(out.status.success());
}
