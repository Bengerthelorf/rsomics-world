use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bbduk"))
}

const ADAPTER: &str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";

// 40 bp clean prefix keeps the survivor ≥ minlength
#[test]
fn kfilter_removes_contaminant_keeps_clean() {
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.fq");
    let clean = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let dirty = format!("ACGTACGTACGTACGTACGT{ADAPTER}");
    let q_clean = "I".repeat(clean.len());
    let q_dirty = "I".repeat(dirty.len());
    fs::write(
        &inp,
        format!("@clean\n{clean}\n+\n{q_clean}\n@dirty\n{dirty}\n+\n{q_dirty}\n"),
    )
    .unwrap();

    let out = Command::new(bin())
        .args(["--literal", ADAPTER, "-i"])
        .arg(&inp)
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    assert!(s.contains("@clean"), "clean read must survive: {s:?}");
    assert!(!s.contains("@dirty"), "contaminant must be removed: {s:?}");
}

#[test]
fn ktrim_right_trims_adapter() {
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.fq");
    let prefix = "ACGTACGTACGTACGTACGT";
    let read = format!("{prefix}{ADAPTER}");
    let q = "I".repeat(read.len());
    fs::write(&inp, format!("@r\n{read}\n+\n{q}\n")).unwrap();

    let out = Command::new(bin())
        .args(["--literal", ADAPTER, "--ktrim", "r", "-i"])
        .arg(&inp)
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    let seq = s.lines().nth(1).unwrap_or("");
    assert_eq!(seq, prefix, "adapter+after must be trimmed: got {seq:?}");
}

#[test]
fn paired_reads_dropped_together() {
    let tmp = tempfile::tempdir().unwrap();
    let r1 = tmp.path().join("r1.fq");
    let r2 = tmp.path().join("r2.fq");
    let o1 = tmp.path().join("o1.fq");
    let o2 = tmp.path().join("o2.fq");
    let clean = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let dirty = format!("ACGTACGTACGTACGTACGT{ADAPTER}");
    fs::write(
        &r1,
        format!("@p1\n{clean}\n+\n{}\n", "I".repeat(clean.len())),
    )
    .unwrap();
    fs::write(
        &r2,
        format!("@p2\n{dirty}\n+\n{}\n", "I".repeat(dirty.len())),
    )
    .unwrap();

    let out = Command::new(bin())
        .args(["--literal", ADAPTER, "-i"])
        .arg(&r1)
        .args(["-I"])
        .arg(&r2)
        .args(["-o"])
        .arg(&o1)
        .args(["-O"])
        .arg(&o2)
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(
        fs::read_to_string(&o1).unwrap().is_empty(),
        "R1 must be dropped because its mate matched"
    );
    assert!(fs::read_to_string(&o2).unwrap().is_empty());
}

#[test]
fn mismatched_pair_counts_fail_loud() {
    let tmp = tempfile::tempdir().unwrap();
    let r1 = tmp.path().join("r1.fq");
    let r2 = tmp.path().join("r2.fq");
    let o1 = tmp.path().join("o1.fq");
    let o2 = tmp.path().join("o2.fq");
    fs::write(&r1, "@a\nACGT\n+\nIIII\n@b\nACGT\n+\nIIII\n").unwrap();
    let mut f2 = fs::File::create(&r2).unwrap();
    f2.write_all(b"@a\nACGT\n+\nIIII\n").unwrap();
    let out = Command::new(bin())
        .args(["-i"])
        .arg(&r1)
        .args(["-I"])
        .arg(&r2)
        .args(["-o"])
        .arg(&o1)
        .args(["-O"])
        .arg(&o2)
        .output()
        .expect("spawn");
    assert!(!out.status.success(), "unequal R1/R2 counts must fail loud");
}

#[test]
fn invalid_k_fails_loud() {
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.fq");
    fs::write(&inp, "@r\nACGT\n+\nIIII\n").unwrap();
    let out = Command::new(bin())
        .args(["--k", "40", "--literal", "ACGT", "-i"])
        .arg(&inp)
        .output()
        .expect("spawn");
    assert!(
        !out.status.success(),
        "--k 40 exceeds the 2-bit codec limit (31) and must fail loud"
    );
}
