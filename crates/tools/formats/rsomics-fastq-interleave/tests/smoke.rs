use std::path::PathBuf;
use std::process::Command;
fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-interleave"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}
#[test]
fn interleaves_pairs() {
    let out = Command::new(bin())
        .args(["-i"])
        .arg(fixture("r1.fq"))
        .args(["-I"])
        .arg(fixture("r2.fq"))
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    let names: Vec<&str> = s.lines().filter(|l| l.starts_with('@')).collect();
    assert_eq!(
        names,
        vec!["@r1/1", "@r1/2", "@r2/1", "@r2/2"],
        "R1/R2 interleaved: {names:?}"
    );
}
