use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-merge"))
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

fn bedtools_merge(input: &std::path::Path) -> Vec<u8> {
    // bedtools merge requires sorted input.
    let sorted = Command::new("bedtools")
        .args(["sort", "-i", input.to_str().unwrap()])
        .output()
        .expect("bedtools sort");
    assert!(
        sorted.status.success(),
        "{}",
        String::from_utf8_lossy(&sorted.stderr)
    );
    let mut child = std::process::Command::new("bedtools")
        .args(["merge", "-i", "-"])
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("bedtools merge spawn");
    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(&sorted.stdout)
        .unwrap();
    let out = child.wait_with_output().expect("bedtools merge wait");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

#[test]
fn merge_matches_bedtools_on_overlap_set() {
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
            "chr1\t100\t200",
            "chr1\t150\t250",
            "chr1\t300\t400",
            "chr2\t100\t200",
            "chr2\t300\t400",
            "chr2\t350\t450",
        ],
    );
    let status = Command::new(ours())
        .args(["-i", inp.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .status()
        .unwrap();
    assert!(status.success());
    let ours_bytes = std::fs::read(&out).unwrap();
    let theirs_bytes = bedtools_merge(&inp);
    assert_eq!(
        ours_bytes,
        theirs_bytes,
        "rsomics: {}\nbedtools: {}",
        String::from_utf8_lossy(&ours_bytes),
        String::from_utf8_lossy(&theirs_bytes)
    );
}

#[test]
fn merge_matches_bedtools_on_touching_intervals() {
    assert!(
        bedtools_available(),
        "compat test requires bedtools on PATH"
    );
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.bed");
    let out = tmp.path().join("ours.bed");
    // [100,200)+[200,300) — bedtools merge -d 0 collapses these.
    write_bed(&inp, &["chr1\t100\t200", "chr1\t200\t300"]);
    let status = Command::new(ours())
        .args(["-i", inp.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .status()
        .unwrap();
    assert!(status.success());
    let ours_bytes = std::fs::read(&out).unwrap();
    let theirs_bytes = bedtools_merge(&inp);
    assert_eq!(ours_bytes, theirs_bytes);
}
