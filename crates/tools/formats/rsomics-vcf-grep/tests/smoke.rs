use std::path::PathBuf;
use std::process::Command;
fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf-grep"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}
#[test]
fn grep_rsid() {
    let out = Command::new(bin())
        .args(["-p", "rs[0-9]+"])
        .arg(fixture("two.vcf"))
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    let count = s
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .count();
    assert_eq!(count, 1, "1 variant has rsID: {s}");
}
