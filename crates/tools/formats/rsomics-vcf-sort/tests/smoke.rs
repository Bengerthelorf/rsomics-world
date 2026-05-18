use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf-sort"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn sorts_by_chrom_pos() {
    let out = Command::new(bin())
        .arg(fixture("unsorted.vcf"))
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    let variants: Vec<&str> = s
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert!(
        variants[0].starts_with("chr1"),
        "chr1 first: {}",
        variants[0]
    );
    assert!(
        variants[1].starts_with("chr2"),
        "chr2 second: {}",
        variants[1]
    );
}
