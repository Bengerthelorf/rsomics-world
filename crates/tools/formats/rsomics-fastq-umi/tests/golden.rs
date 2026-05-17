// Hand-computed fastp-0.20.1 UMI golden fixtures: gates correctness where fastp is not installed.
use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-umi"))
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn se_umi_len6_golden() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("o.fq");
    let st = Command::new(ours())
        .args([
            "-i",
            fixture("se_umi.fastq").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--umi_len",
            "6",
        ])
        .status()
        .unwrap();
    assert!(st.success());

    let expected = "\
@readA:AACCGG 1:N:0:ACGT
TTACGTACGTACGT
+
IIIIIIIIIIIIII
@readB:TTGGCC
AATGCATGCATGCA
+
FFFFFFFFFFFFFF
@readC:GGGGCC desc here
CCAAAATTTTACGT
+
HHHHHHHHHHHHHH
";
    assert_eq!(std::fs::read_to_string(&out).unwrap(), expected);
}

#[test]
fn se_umi_len5_skip2_golden() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("o.fq");
    let st = Command::new(ours())
        .args([
            "-i",
            fixture("se_umi.fastq").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--umi_len",
            "5",
            "--umi_skip",
            "2",
        ])
        .status()
        .unwrap();
    assert!(st.success());

    let expected = "\
@readA:AACCG 1:N:0:ACGT
TACGTACGTACGT
+
IIIIIIIIIIIII
@readB:TTGGC
ATGCATGCATGCA
+
FFFFFFFFFFFFF
@readC:GGGGC desc here
CAAAATTTTACGT
+
HHHHHHHHHHHHH
";
    assert_eq!(std::fs::read_to_string(&out).unwrap(), expected);
}

// fastp 0.20.1 Read::trimFront clamps to length()-1 — reads shorter than umi_len keep their last base.
#[test]
fn se_umi_short_len8_golden() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("o.fq");
    let st = Command::new(ours())
        .args([
            "-i",
            fixture("se_umi_short.fastq").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--umi_len",
            "8",
        ])
        .status()
        .unwrap();
    assert!(st.success());

    let expected = "\
@shortA:AACCG
G
+
I
@shortB:TTG
G
+
F
@shortC:A
A
+
I
";
    assert_eq!(std::fs::read_to_string(&out).unwrap(), expected);
}

// index1: UMI from firstIndex (backward scan to last ':'/'+''), stamped without trimming seq/qual; no index field → pass-through (fastp if(!umi.empty()) guard).
#[test]
fn se_umi_index1_golden() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("o.fq");
    let st = Command::new(ours())
        .args([
            "-i",
            fixture("se_umi.fastq").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--umi_loc",
            "index1",
        ])
        .status()
        .unwrap();
    assert!(st.success());

    let expected = "\
@readA:ACGT 1:N:0:ACGT
AACCGGTTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
@readB
TTGGCCAATGCATGCATGCA
+
FFFFFFFFFFFFFFFFFFFF
@readC desc here
GGGGCCCCAAAATTTTACGT
+
HHHHHHHHHHHHHHHHHHHH
";
    assert_eq!(std::fs::read_to_string(&out).unwrap(), expected);
}

// fastp 0.20.1 trimFront throws on zero-length reads; fail loud rather than fabricate a record.
#[test]
fn empty_source_read_errors_cli() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("o.fq");
    let st = Command::new(ours())
        .args([
            "-i",
            fixture("se_umi_empty.fastq").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--umi_len",
            "4",
        ])
        .status()
        .unwrap();
    assert!(
        !st.success(),
        "empty UMI source read must fail loud, not emit a fabricated record"
    );
}

// fastp 0.20.1 rejects PE-only locations (read2/index2/per_index/per_read) without paired input at validation time.
#[test]
fn pe_only_loc_in_se_rejected_cli() {
    for loc in ["per_index", "per_read", "read2", "index2"] {
        let tmp = tempfile::tempdir().unwrap();
        let out = tmp.path().join("o.fq");
        let st = Command::new(ours())
            .args([
                "-i",
                fixture("se_umi.fastq").to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                "--umi_loc",
                loc,
                "--umi_len",
                "4",
            ])
            .status()
            .unwrap();
        assert!(
            !st.success(),
            "--umi_loc {loc} in SE must fail loud (fastp rejects it), got success"
        );
    }
}
