use std::path::PathBuf;
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bam-idxstats"))
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

#[test]
fn idxstats_matches_samtools() {
    if !samtools_available() {
        eprintln!("samtools not on PATH — skipping compat test");
        return;
    }
    let input = fixture("small.bam");
    if !input.exists() {
        eprintln!("golden fixture small.bam missing — skipping");
        return;
    }

    let ours_out = Command::new(ours())
        .arg(&input)
        .output()
        .expect("spawn ours");
    assert!(
        ours_out.status.success(),
        "rsomics-bam-idxstats failed: {}",
        String::from_utf8_lossy(&ours_out.stderr)
    );

    let theirs_out = Command::new("samtools")
        .args(["idxstats", input.to_str().unwrap()])
        .output()
        .expect("spawn samtools");
    assert!(
        theirs_out.status.success(),
        "samtools idxstats failed: {}",
        String::from_utf8_lossy(&theirs_out.stderr)
    );

    let ours = String::from_utf8(ours_out.stdout).unwrap();
    let theirs = String::from_utf8(theirs_out.stdout).unwrap();

    assert_eq!(
        ours.trim(),
        theirs.trim(),
        "output differs:\n--- ours ---\n{ours}\n--- samtools ---\n{theirs}"
    );
}
