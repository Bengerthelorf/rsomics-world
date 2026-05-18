use std::path::PathBuf;
use std::process::Command;
fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-sample"))
}
#[test]
fn help_works() {
    let out = Command::new(ours())
        .args(["--help"])
        .output()
        .expect("spawn");
    assert!(out.status.success());
}
