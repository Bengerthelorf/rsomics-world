use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-sort"))
}

fn write_bed(path: &std::path::Path, lines: &[&str]) {
    let mut f = std::fs::File::create(path).unwrap();
    for line in lines {
        writeln!(f, "{line}").unwrap();
    }
}

#[test]
fn sorts_within_chrom_by_start() {
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.bed");
    let out = tmp.path().join("out.bed");
    write_bed(
        &inp,
        &["chr1\t300\t400", "chr1\t100\t200", "chr1\t200\t300"],
    );
    let status = Command::new(ours())
        .args(["-i", inp.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(out).unwrap();
    assert_eq!(got, "chr1\t100\t200\nchr1\t200\t300\nchr1\t300\t400\n");
}

#[test]
fn lexicographic_chrom_order_matches_bedtools() {
    // bedtools' default sort is lexicographic, so "chr10" sorts before "chr2".
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.bed");
    let out = tmp.path().join("out.bed");
    write_bed(
        &inp,
        &["chr2\t100\t200", "chr10\t100\t200", "chr1\t100\t200"],
    );
    let status = Command::new(ours())
        .args(["-i", inp.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(out).unwrap();
    assert_eq!(got, "chr1\t100\t200\nchr10\t100\t200\nchr2\t100\t200\n");
}
