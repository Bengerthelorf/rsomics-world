use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-gff-filter"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn filter_by_type_exon() {
    let out = Command::new(bin())
        .args(["--type", "exon"])
        .arg(fixture("small.gff"))
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    let feature_count = s
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .count();
    assert_eq!(feature_count, 2, "2 exons: {s}");
}

#[test]
fn filter_by_type_gene() {
    let out = Command::new(bin())
        .args(["--type", "gene"])
        .arg(fixture("small.gff"))
        .output()
        .expect("spawn");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let count = s
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .count();
    assert_eq!(count, 2, "2 genes: {s}");
}
