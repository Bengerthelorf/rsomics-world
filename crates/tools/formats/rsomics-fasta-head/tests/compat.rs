use std::path::PathBuf;
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-head"))
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
fn head_5_matches_seqkit() {
    if !seqkit_available() {
        eprintln!("seqkit not on PATH — skipping compat test");
        return;
    }
    let input = fixture("twenty.fa");

    let ours_out = Command::new(ours())
        .args(["-n", "5"])
        .arg(&input)
        .output()
        .expect("spawn ours");
    assert!(ours_out.status.success());

    let seqkit_out = Command::new("seqkit")
        .args(["head", "-n", "5"])
        .arg(&input)
        .output()
        .expect("spawn seqkit");
    assert!(seqkit_out.status.success());

    let ours_names: Vec<String> = String::from_utf8(ours_out.stdout)
        .unwrap()
        .lines()
        .filter(|l| l.starts_with('>'))
        .map(String::from)
        .collect();
    let seqkit_names: Vec<String> = String::from_utf8(seqkit_out.stdout)
        .unwrap()
        .lines()
        .filter(|l| l.starts_with('>'))
        .map(String::from)
        .collect();

    assert_eq!(ours_names.len(), 5);
    assert_eq!(ours_names, seqkit_names, "names must match seqkit head");
}
