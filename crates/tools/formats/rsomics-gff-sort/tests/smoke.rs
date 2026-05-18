use std::path::PathBuf;
use std::process::Command;
fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-gff-sort"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}
#[test]
fn sorts_by_chrom_pos() {
    let out = Command::new(bin())
        .arg(fixture("unsorted.gff"))
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    let first_data = s
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .next()
        .unwrap();
    assert!(first_data.starts_with("chr1"), "chr1 first: {first_data}");
}
