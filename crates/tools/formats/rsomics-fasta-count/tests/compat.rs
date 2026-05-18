use std::path::PathBuf;
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-count"))
}

fn seqkit_available() -> bool {
    Command::new("seqkit")
        .arg("version")
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
fn count_matches_seqkit() {
    if !seqkit_available() {
        eprintln!("seqkit not on PATH — skipping compat test");
        return;
    }
    let input = fixture("five.fa");

    let ours_out = Command::new(ours())
        .arg(&input)
        .output()
        .expect("spawn ours");
    assert!(ours_out.status.success());
    let our_count: u64 = String::from_utf8(ours_out.stdout)
        .unwrap()
        .trim()
        .parse()
        .unwrap();

    let seqkit_out = Command::new("seqkit")
        .args(["stats", "-T"])
        .arg(&input)
        .output()
        .expect("spawn seqkit");
    assert!(seqkit_out.status.success());
    let seqkit_str = String::from_utf8(seqkit_out.stdout).unwrap();
    let seqkit_count: u64 = seqkit_str
        .lines()
        .nth(1)
        .and_then(|l| l.split('\t').nth(3))
        .and_then(|s| s.replace(',', "").parse().ok())
        .unwrap_or(0);

    assert_eq!(
        our_count, seqkit_count,
        "count mismatch: ours={our_count} seqkit={seqkit_count}"
    );
}
