use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-merge"))
}

// fastp v0.20.1's overlapanalysis test() pair — known to overlap
// (offset 10, overlap_len 79), so the binary must emit one merged read.
#[test]
fn emits_a_merged_read() {
    let tmp = tempfile::tempdir().unwrap();
    let r1 = tmp.path().join("r1.fq");
    let r2 = tmp.path().join("r2.fq");
    let q: String = "I".repeat(89);
    fs::write(
        &r1,
        format!("@p/1\nCAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGC\n+\n{q}\n"),
    )
    .unwrap();
    fs::write(
        &r2,
        format!("@p/2\nACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGTCC\n+\n{q}\n"),
    )
    .unwrap();

    let out = Command::new(bin())
        .args(["--in1"])
        .arg(&r1)
        .args(["--in2"])
        .arg(&r2)
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    assert!(s.contains(" merged_"), "expected a merged read header: {s}");
    assert_eq!(s.lines().count(), 4, "exactly one merged FASTQ record");
}

#[test]
fn mismatched_pair_counts_fail_loud() {
    let tmp = tempfile::tempdir().unwrap();
    let r1 = tmp.path().join("r1.fq");
    let r2 = tmp.path().join("r2.fq");
    fs::write(&r1, "@a\nACGT\n+\nIIII\n@b\nACGT\n+\nIIII\n").unwrap();
    let mut f2 = fs::File::create(&r2).unwrap();
    f2.write_all(b"@a\nACGT\n+\nIIII\n").unwrap();
    let out = Command::new(bin())
        .args(["--in1"])
        .arg(&r1)
        .args(["--in2"])
        .arg(&r2)
        .output()
        .expect("spawn");
    assert!(!out.status.success(), "unequal R1/R2 counts must fail loud");
}
