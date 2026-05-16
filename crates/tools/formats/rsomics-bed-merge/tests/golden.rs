use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-merge"))
}

fn write_bed(path: &std::path::Path, lines: &[&str]) {
    let mut f = std::fs::File::create(path).unwrap();
    for line in lines {
        writeln!(f, "{line}").unwrap();
    }
}

#[test]
fn merges_two_overlapping_intervals() {
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.bed");
    let out = tmp.path().join("out.bed");
    write_bed(&inp, &["chr1\t100\t200", "chr1\t150\t250"]);
    let status = Command::new(ours())
        .args(["-i", inp.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(out).unwrap();
    assert_eq!(got, "chr1\t100\t250\n");
}

#[test]
fn keeps_non_overlapping_intervals_separate() {
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.bed");
    let out = tmp.path().join("out.bed");
    write_bed(&inp, &["chr1\t100\t200", "chr1\t300\t400"]);
    let status = Command::new(ours())
        .args(["-i", inp.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(out).unwrap();
    assert_eq!(got, "chr1\t100\t200\nchr1\t300\t400\n");
}

#[test]
fn merges_touching_intervals() {
    // [100,200) + [200,300) — bedtools merge -d 0 collapses to [100,300).
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.bed");
    let out = tmp.path().join("out.bed");
    write_bed(&inp, &["chr1\t100\t200", "chr1\t200\t300"]);
    let status = Command::new(ours())
        .args(["-i", inp.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(out).unwrap();
    assert_eq!(got, "chr1\t100\t300\n");
}

#[test]
fn multi_chrom_independent_merge() {
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.bed");
    let out = tmp.path().join("out.bed");
    write_bed(
        &inp,
        &[
            "chr1\t100\t200",
            "chr1\t150\t250",
            "chr2\t100\t200",
            "chr2\t300\t400",
        ],
    );
    let status = Command::new(ours())
        .args(["-i", inp.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(out).unwrap();
    assert_eq!(got, "chr1\t100\t250\nchr2\t100\t200\nchr2\t300\t400\n");
}

#[test]
fn unsorted_input_fails_loud() {
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.bed");
    let out = tmp.path().join("out.bed");
    write_bed(&inp, &["chr1\t300\t400", "chr1\t100\t200"]);
    let result = Command::new(ours())
        .args(["-i", inp.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(
        String::from_utf8_lossy(&result.stderr).contains("not sorted"),
        "stderr: {}",
        String::from_utf8_lossy(&result.stderr)
    );
}

#[test]
fn reappearing_chrom_fails_loud() {
    // chr1 -> chr2 -> chr1: bedtools errors ("out of order"); we must too,
    // not silently emit an interleaved split that downstream set-algebra
    // would mis-handle.
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.bed");
    let out = tmp.path().join("out.bed");
    write_bed(
        &inp,
        &["chr1\t100\t200", "chr2\t100\t200", "chr1\t300\t400"],
    );
    let result = Command::new(ours())
        .args(["-i", inp.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(
        String::from_utf8_lossy(&result.stderr).contains("not sorted"),
        "stderr: {}",
        String::from_utf8_lossy(&result.stderr)
    );
}
