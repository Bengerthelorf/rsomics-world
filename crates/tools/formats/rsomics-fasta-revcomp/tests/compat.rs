use std::path::PathBuf;
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-revcomp"))
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
fn revcomp_matches_seqkit() {
    if !seqkit_available() {
        eprintln!("seqkit not on PATH — skipping compat test");
        return;
    }
    let input = fixture("one.fa");

    let ours_out = Command::new(ours())
        .arg(&input)
        .output()
        .expect("spawn ours");
    assert!(ours_out.status.success());

    let seqkit_out = Command::new("seqkit")
        .args(["seq", "-r", "-p"])
        .arg(&input)
        .output()
        .expect("spawn seqkit");
    assert!(seqkit_out.status.success());

    let ours_seq = String::from_utf8(ours_out.stdout)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('>'))
        .collect::<Vec<_>>()
        .join("");
    let seqkit_seq = String::from_utf8(seqkit_out.stdout)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('>'))
        .collect::<Vec<_>>()
        .join("");

    assert_eq!(
        ours_seq, seqkit_seq,
        "revcomp sequences must match: ours={ours_seq} seqkit={seqkit_seq}"
    );
}
