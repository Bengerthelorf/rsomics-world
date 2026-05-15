use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-intersect"))
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

fn bedtools_intersect(a: &std::path::Path, b: &std::path::Path) -> Vec<u8> {
    let out = Command::new("bedtools")
        .args([
            "intersect",
            "-a",
            a.to_str().unwrap(),
            "-b",
            b.to_str().unwrap(),
        ])
        .output()
        .expect("bedtools intersect");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

#[test]
fn intersect_matches_bedtools() {
    assert!(
        bedtools_available(),
        "compat test requires bedtools on PATH"
    );
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("a.bed");
    let b = tmp.path().join("b.bed");
    let out = tmp.path().join("ours.bed");
    write_bed(&a, &["chr1\t100\t200", "chr1\t300\t400", "chr2\t100\t300"]);
    write_bed(&b, &["chr1\t150\t350", "chr2\t200\t250"]);
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
    // bedtools' default intersect preserves input order from `a`; ours
    // chrom-sorts on output. Normalise both by sorting before compare.
    let mut ours_lines: Vec<String> = std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .map(String::from)
        .collect();
    let mut theirs_lines: Vec<String> = String::from_utf8(bedtools_intersect(&a, &b))
        .unwrap()
        .lines()
        .map(String::from)
        .collect();
    ours_lines.sort();
    theirs_lines.sort();
    assert_eq!(ours_lines, theirs_lines);
}
