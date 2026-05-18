use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf-stats"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn counts_snps_and_indels() {
    let out = Command::new(bin())
        .arg(fixture("small.vcf"))
        .output()
        .expect("spawn");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    assert!(s.contains("Total variants:\t4"), "expected 4 total: {s}");
    assert!(
        s.contains("SNPs:\t3"),
        "expected 3 SNPs (A>G, C>T, A>C): {s}"
    );
    assert!(
        s.contains("Deletions:\t1"),
        "expected 1 deletion (AT>A): {s}"
    );
}

#[test]
fn titv_ratio_correct() {
    let out = Command::new(bin())
        .arg(fixture("small.vcf"))
        .output()
        .expect("spawn");
    let s = String::from_utf8(out.stdout).unwrap();
    // A>G (ti), C>T (ti), A>C (tv) → Ti=2, Tv=1, Ti/Tv=2.0
    assert!(s.contains("Transitions:\t2"), "expected 2 transitions: {s}");
    assert!(
        s.contains("Transversions:\t1"),
        "expected 1 transversion: {s}"
    );
    assert!(
        s.contains("Ti/Tv ratio:\t2.0000"),
        "expected Ti/Tv=2.0: {s}"
    );
}
