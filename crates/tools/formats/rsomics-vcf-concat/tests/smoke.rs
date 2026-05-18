use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf-concat"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn concatenates_two_vcfs() {
    let out = Command::new(bin())
        .arg(fixture("a.vcf"))
        .arg(fixture("b.vcf"))
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    let variant_count = s
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .count();
    assert_eq!(variant_count, 2, "should have both variants: {s}");
    assert!(s.contains("chr1\t100"), "has chr1 variant");
    assert!(s.contains("chr2\t200"), "has chr2 variant");
}
