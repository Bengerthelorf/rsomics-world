use std::path::PathBuf;
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bam-flagstat"))
}

fn samtools_available() -> bool {
    Command::new("samtools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn run_ours(input: &std::path::Path) -> String {
    let out = Command::new(ours())
        .arg(input)
        .output()
        .expect("spawn ours");
    assert!(
        out.status.success(),
        "rsomics-bam-flagstat failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn run_samtools(input: &std::path::Path) -> String {
    let out = Command::new("samtools")
        .args(["flagstat", input.to_str().unwrap()])
        .output()
        .expect("spawn samtools");
    assert!(
        out.status.success(),
        "samtools flagstat failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

#[test]
fn flagstat_matches_samtools() {
    if !samtools_available() {
        eprintln!("samtools not on PATH — skipping compat test");
        return;
    }
    let input = fixture("small.bam");
    if !input.exists() {
        eprintln!("golden fixture small.bam missing — skipping compat test");
        return;
    }

    let ours = run_ours(&input);
    let theirs = run_samtools(&input);

    assert_eq!(
        ours.trim(),
        theirs.trim(),
        "output differs:\n--- ours ---\n{ours}\n--- samtools ---\n{theirs}"
    );
}
