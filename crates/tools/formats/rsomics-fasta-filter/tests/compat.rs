use std::path::PathBuf;
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-filter"))
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
fn min_len_matches_seqkit() {
    if !seqkit_available() {
        eprintln!("seqkit not on PATH — skipping compat test");
        return;
    }
    let input = fixture("mixed.fa");

    let ours_out = Command::new(ours())
        .args(["-m", "20"])
        .arg(&input)
        .output()
        .expect("spawn ours");
    assert!(ours_out.status.success());

    let seqkit_out = Command::new("seqkit")
        .args(["seq", "-m", "20"])
        .arg(&input)
        .output()
        .expect("spawn seqkit");
    assert!(seqkit_out.status.success());

    let ours_count = String::from_utf8(ours_out.stdout)
        .unwrap()
        .lines()
        .filter(|l| l.starts_with('>'))
        .count();
    let seqkit_count = String::from_utf8(seqkit_out.stdout)
        .unwrap()
        .lines()
        .filter(|l| l.starts_with('>'))
        .count();

    assert_eq!(
        ours_count, seqkit_count,
        "min-len filter count: ours={ours_count} seqkit={seqkit_count}"
    );
}
