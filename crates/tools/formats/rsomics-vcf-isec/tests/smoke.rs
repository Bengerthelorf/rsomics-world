use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf-isec"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn shared_variant() {
    let out = Command::new(bin())
        .args(["-a"])
        .arg(fixture("a.vcf"))
        .args(["-b"])
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
    assert_eq!(variant_count, 1, "chr1:100 A>G is shared: {s}");
}
