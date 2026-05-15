use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-subtract"))
}

fn write_bed(path: &std::path::Path, lines: &[&str]) {
    let mut f = std::fs::File::create(path).unwrap();
    for line in lines {
        writeln!(f, "{line}").unwrap();
    }
}

#[test]
fn punches_holes_in_a() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("a.bed");
    let b = tmp.path().join("b.bed");
    let out = tmp.path().join("out.bed");
    write_bed(&a, &["chr1\t100\t500"]);
    write_bed(&b, &["chr1\t200\t300", "chr1\t400\t450"]);
    let status = Command::new(ours())
        .args([
            "-a",
            a.to_str().unwrap(),
            "-b",
            b.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(out).unwrap();
    assert_eq!(got, "chr1\t100\t200\nchr1\t300\t400\nchr1\t450\t500\n");
}

#[test]
fn full_cover_yields_empty() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("a.bed");
    let b = tmp.path().join("b.bed");
    let out = tmp.path().join("out.bed");
    write_bed(&a, &["chr1\t100\t200"]);
    write_bed(&b, &["chr1\t50\t250"]);
    let status = Command::new(ours())
        .args([
            "-a",
            a.to_str().unwrap(),
            "-b",
            b.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(out).unwrap();
    assert_eq!(got, "");
}

#[test]
fn b_in_different_chrom_is_noop() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("a.bed");
    let b = tmp.path().join("b.bed");
    let out = tmp.path().join("out.bed");
    write_bed(&a, &["chr1\t100\t200"]);
    write_bed(&b, &["chr2\t100\t200"]);
    let status = Command::new(ours())
        .args([
            "-a",
            a.to_str().unwrap(),
            "-b",
            b.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(out).unwrap();
    assert_eq!(got, "chr1\t100\t200\n");
}
