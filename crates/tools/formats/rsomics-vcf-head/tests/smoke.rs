use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf-head"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn head_3_variants() {
    let out = Command::new(bin())
        .args(["-n", "3"])
        .arg(fixture("five.vcf"))
        .output()
        .expect("spawn");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let variant_count = s
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .count();
    assert_eq!(variant_count, 3, "expected 3 variants: {s}");
}

#[test]
fn preserves_header() {
    let out = Command::new(bin())
        .args(["-n", "1"])
        .arg(fixture("five.vcf"))
        .output()
        .expect("spawn");
    let s = String::from_utf8(out.stdout).unwrap();
    assert!(s.contains("##fileformat"), "header must be preserved");
    assert!(s.contains("#CHROM"), "column header must be preserved");
}
