use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-complement"))
}

fn write_lines(path: &std::path::Path, lines: &[&str]) {
    let mut f = std::fs::File::create(path).unwrap();
    for line in lines {
        writeln!(f, "{line}").unwrap();
    }
}

#[test]
fn complement_yields_inter_interval_gaps() {
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.bed");
    let g = tmp.path().join("genome.tsv");
    let out = tmp.path().join("out.bed");
    write_lines(&inp, &["chr1\t100\t200", "chr1\t400\t500"]);
    write_lines(&g, &["chr1\t1000"]);
    let status = Command::new(ours())
        .args([
            "-i",
            inp.to_str().unwrap(),
            "-g",
            g.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(out).unwrap();
    assert_eq!(got, "chr1\t0\t100\nchr1\t200\t400\nchr1\t500\t1000\n");
}

#[test]
fn empty_input_yields_full_genome() {
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.bed");
    let g = tmp.path().join("genome.tsv");
    let out = tmp.path().join("out.bed");
    write_lines(&inp, &[]);
    write_lines(&g, &["chr1\t1000", "chr2\t500"]);
    let status = Command::new(ours())
        .args([
            "-i",
            inp.to_str().unwrap(),
            "-g",
            g.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(out).unwrap();
    assert_eq!(got, "chr1\t0\t1000\nchr2\t0\t500\n");
}

#[test]
fn full_coverage_yields_no_complement() {
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.bed");
    let g = tmp.path().join("genome.tsv");
    let out = tmp.path().join("out.bed");
    write_lines(&inp, &["chr1\t0\t1000"]);
    write_lines(&g, &["chr1\t1000"]);
    let status = Command::new(ours())
        .args([
            "-i",
            inp.to_str().unwrap(),
            "-g",
            g.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(out).unwrap();
    assert_eq!(got, "");
}
