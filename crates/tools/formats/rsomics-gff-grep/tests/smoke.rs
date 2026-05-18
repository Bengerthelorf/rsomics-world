use std::path::PathBuf;
use std::process::Command;
fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-gff-grep"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}
#[test]
fn grep_tp53() {
    let out = Command::new(bin())
        .args(["-p", "TP53"])
        .arg(fixture("two.gff"))
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
    assert_eq!(count, 1, "1 TP53 feature: {s}");
}
