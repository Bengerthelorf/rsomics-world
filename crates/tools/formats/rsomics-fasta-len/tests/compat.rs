use std::path::PathBuf;
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-len"))
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
fn lengths_match_seqkit() {
    if !seqkit_available() {
        eprintln!("seqkit not on PATH — skipping compat test");
        return;
    }
    let input = fixture("three.fa");

    let ours_out = Command::new(ours())
        .args(["--tab"])
        .arg(&input)
        .output()
        .expect("spawn ours");
    assert!(ours_out.status.success());

    let seqkit_out = Command::new("seqkit")
        .args(["fx2tab", "-n", "-l"])
        .arg(&input)
        .output()
        .expect("spawn seqkit");
    assert!(seqkit_out.status.success());

    let ours_lines: Vec<String> = String::from_utf8(ours_out.stdout)
        .unwrap()
        .lines()
        .map(|l| l.trim().to_string())
        .collect();
    let seqkit_lines: Vec<String> = String::from_utf8(seqkit_out.stdout)
        .unwrap()
        .lines()
        .map(|l| l.trim().to_string())
        .collect();

    assert_eq!(
        ours_lines.len(),
        seqkit_lines.len(),
        "record count mismatch"
    );
    for (o, s) in ours_lines.iter().zip(seqkit_lines.iter()) {
        let o_parts: Vec<&str> = o.split('\t').collect();
        let s_parts: Vec<&str> = s.split('\t').collect();
        assert_eq!(
            o_parts.last(),
            s_parts.last(),
            "length mismatch: ours={o} seqkit={s}"
        );
    }
}
