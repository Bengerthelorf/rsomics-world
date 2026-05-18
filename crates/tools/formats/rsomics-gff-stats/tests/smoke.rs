use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-gff-stats"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn counts_features() {
    let out = Command::new(bin())
        .arg(fixture("small.gff"))
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    assert!(s.contains("Total features:\t3"), "3 features: {s}");
    assert!(s.contains("gene\t2"), "2 genes: {s}");
    assert!(s.contains("exon\t1"), "1 exon: {s}");
}
