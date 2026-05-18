use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-gc"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn gc_content_correct() {
    let out = Command::new(bin())
        .arg(fixture("two.fa"))
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    assert!(s.contains("100.00"), "GCGCGCGC = 100% GC: {s}");
    assert!(s.contains("0.00"), "AAAAAAAA = 0% GC: {s}");
}
