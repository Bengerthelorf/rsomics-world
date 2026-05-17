// BFC byte output depends on its counting-hash collision profile → byte diff is undefined
// across counter impls. Test: on a high-coverage golden fixture, ours and bfc must agree on
// the corrected base (implementation-independent). Skipped when bfc not on PATH; authoritative on 4090/CI.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-correct"))
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn bfc_available() -> bool {
    // exit status varies by bfc build; "binary ran at all" is the signal.
    Command::new("bfc")
        .arg("-v")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok()
}

fn seqs(bytes: &[u8]) -> Vec<Vec<u8>> {
    let text = String::from_utf8_lossy(bytes);
    text.lines()
        .collect::<Vec<_>>()
        .chunks(4)
        .filter(|c| c.len() == 4)
        .map(|c| c[1].to_ascii_uppercase().into_bytes())
        .collect()
}

#[test]
fn single_substitution_outcome_matches_bfc() {
    if !bfc_available() {
        eprintln!(
            "SKIP: bfc not on PATH — semantic compat oracle unavailable (authoritative on 4090/CI)"
        );
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let ours_out = tmp.path().join("ours.fq");
    let theirs_out = tmp.path().join("theirs.fq");
    let input = fixture("one_subst.fastq");

    let st = Command::new(ours())
        .args([
            "-i",
            input.to_str().unwrap(),
            "-o",
            ours_out.to_str().unwrap(),
            "-k",
            "17",
            "-c",
            "3",
        ])
        .status()
        .unwrap();
    assert!(st.success(), "ours exited non-zero");

    let bfc = Command::new("bfc")
        .args(["-s", "1k", "-k", "17", "-t", "1", input.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(bfc.status.success(), "bfc exited non-zero");
    std::fs::write(&theirs_out, &bfc.stdout).unwrap();

    let a = seqs(&std::fs::read(&ours_out).unwrap());
    let b = seqs(&bfc.stdout);
    assert_eq!(a.len(), b.len(), "record count diverged from bfc");
    // The planted error is at the last record; in a high-coverage context
    // the corrected base is implementation-independent — both must agree
    // on every record's final sequence.
    for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
        assert_eq!(
            x,
            y,
            "record {i}: correction outcome diverged from bfc\n ours: {}\n bfc:  {}",
            String::from_utf8_lossy(x),
            String::from_utf8_lossy(y)
        );
    }
    let _ = Path::new(&theirs_out);
}
