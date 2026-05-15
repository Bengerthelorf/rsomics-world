use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-complement"))
}

fn bedtools_available() -> bool {
    Command::new("bedtools")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
}

fn write_lines(path: &std::path::Path, lines: &[&str]) {
    let mut f = std::fs::File::create(path).unwrap();
    for line in lines {
        writeln!(f, "{line}").unwrap();
    }
}

#[test]
fn complement_matches_bedtools() {
    assert!(
        bedtools_available(),
        "compat test requires bedtools on PATH"
    );
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.bed");
    let g = tmp.path().join("genome.tsv");
    let out = tmp.path().join("ours.bed");
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
    let ours_bytes = std::fs::read(&out).unwrap();
    let sorted = Command::new("bedtools")
        .args(["sort", "-i", inp.to_str().unwrap()])
        .output()
        .expect("bedtools sort");
    assert!(sorted.status.success());
    let mut child = std::process::Command::new("bedtools")
        .args(["complement", "-i", "-", "-g", g.to_str().unwrap()])
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("bedtools complement spawn");
    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(&sorted.stdout)
        .unwrap();
    let theirs = child.wait_with_output().expect("bedtools complement wait");
    assert!(
        theirs.status.success(),
        "{}",
        String::from_utf8_lossy(&theirs.stderr)
    );
    assert_eq!(ours_bytes, theirs.stdout);
}
