//! fastp-independent golden: hand-computed fastp-0.20.1 UMI output so
//! correctness is gated everywhere, not only where fastp is installed.

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

    // readA: umi=AACCGG, name has a space → tag before it; seq/qual lose 6 at 5'
    // readB: umi=TTGGCC, no space → appended
    // readC: umi=GGGGCC, space after "readC" → tag before it
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

    // umi = first 5 bases; trimFront = 5 + 2 = 7 from seq+qual (20 → 13).
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

/// Reads SHORTER than `umi_len` + skip: fastp 0.20.1 `Read::trimFront` clamps
/// to `length()-1`, so every read keeps at least its last base (never emptied).
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

    // shortA len5 umi=AACCG trim min(5,4)=4 → G/I; shortB len3 umi=TTG trim
    // min(3,2)=2 → G/F; shortC len1 umi=A trim min(1,0)=0 → A/I kept.
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

/// `--umi_loc index1`: the UMI is the read-name trailing index field
/// (fastp 0.20.1 `Read::firstIndex`, backward scan to the last `:`/`+`),
/// stamped without trimming seq/qual. A name with no index field is a
/// pass-through (fastp's `if(!umi.empty())` guard), not an error.
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

    // readA name=`readA 1:N:0:ACGT` → firstIndex=ACGT, tag before first space,
    //   seq/qual untouched. readB=`readB` / readC=`readC desc here` have no
    //   index field → empty UMI → pass-through unchanged.
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

/// A zero-length UMI source read has no defined fastp 0.20.1 output (fastp's
/// `trimFront` throws). We fail loud rather than fabricate a record.
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

/// fastp 0.20.1 rejects `read2` / `index2` / `per_index` / `per_read`
/// without a paired input at option-validation time. Match it: SE with a
/// PE-only location must fail loud, not silently emit a partial UMI.
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
