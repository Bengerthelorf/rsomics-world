use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-sort"))
}

fn bedtools_available() -> bool {
    Command::new("bedtools")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
}

fn write_bed(path: &std::path::Path, lines: &[&str]) {
    let mut f = std::fs::File::create(path).unwrap();
    for line in lines {
        writeln!(f, "{line}").unwrap();
    }
}

#[test]
fn sort_matches_bedtools_default() {
    assert!(
        bedtools_available(),
        "compat test requires bedtools on PATH"
    );
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.bed");
    let out = tmp.path().join("ours.bed");
    write_bed(
        &inp,
        &[
            "chr2\t300\t400",
            "chr1\t100\t200",
            "chr10\t50\t60",
            "chr1\t50\t150",
        ],
    );
    let status = Command::new(ours())
        .args(["-i", inp.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .status()
        .unwrap();
    assert!(status.success());
    let ours_bytes = std::fs::read(&out).unwrap();
    let theirs = Command::new("bedtools")
        .args(["sort", "-i", inp.to_str().unwrap()])
        .output()
        .expect("bedtools sort");
    assert!(theirs.status.success());
    assert_eq!(ours_bytes, theirs.stdout);
}
