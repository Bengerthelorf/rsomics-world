use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-filter"))
}

fn run_filter(args: &[&str]) -> Vec<u8> {
    let out = Command::new(ours())
        .args(args)
        .output()
        .expect("spawn rsomics-fastq-filter");
    assert!(
        out.status.success(),
        "rsomics-fastq-filter {:?} failed:\nstdout: {}\nstderr: {}",
        args,
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr),
    );
    out.stdout
}

fn read_fastq_ids(path: &Path) -> Vec<String> {
    let content = std::fs::read_to_string(path).expect("read output");
    content
        .lines()
        .filter(|l| l.starts_with('@'))
        .map(|l| l[1..].split_whitespace().next().unwrap_or("").to_owned())
        .collect()
}

// Quality filter semantics (from fastp src/filter.cpp + src/options.h):
//   unqualified = qual_byte < (qualified_quality_phred + 33)
//   fail iff low_qual_count > unqualified_percent_limit * len / 100.0
//   fail iff n_count > n_base_limit
//
// With defaults: qual_threshold=15 (ASCII 48='0'), unqualified_percent_limit=40, n_base_limit=5
// Fixture reads are 20 bp; limit = 40 * 20 / 100.0 = 8.0; fail iff low_qual_count > 8.
// '!' = ASCII 33 < 48 → unqualified. 'I' = ASCII 73 ≥ 48 → qualified.
#[test]
fn se_quality_filter_pass_fail() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("out.fq");
    let input = fixture("se_quality.fastq");

    run_filter(&[
        "-i",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "-L",
    ]);

    let ids = read_fastq_ids(&out);
    // r1: 0 low-qual → PASS
    assert!(
        ids.contains(&"r1_all_high_qual".to_owned()),
        "r1 should pass"
    );
    // r2: 10 low-qual → 10 > 8.0 → FAIL
    assert!(
        !ids.contains(&"r2_too_many_low_qual".to_owned()),
        "r2 should fail"
    );
    // r3: 8 low-qual → 8 > 8.0 is false → PASS (exact boundary)
    assert!(
        ids.contains(&"r3_exactly_at_limit".to_owned()),
        "r3 at boundary should pass"
    );
    // r4: 9 low-qual → 9 > 8.0 → FAIL (one over boundary)
    assert!(
        !ids.contains(&"r4_one_over_limit".to_owned()),
        "r4 one over limit should fail"
    );
    // r5: 20 low-qual → FAIL
    assert!(
        !ids.contains(&"r5_all_low_qual".to_owned()),
        "r5 should fail"
    );
    // r6: 1 N → 1 > 5 is false → PASS
    assert!(
        ids.contains(&"r6_n_bases_ok".to_owned()),
        "r6 with 1 N should pass"
    );
    // r7: 6 N → 6 > 5 → FAIL
    assert!(
        !ids.contains(&"r7_too_many_n".to_owned()),
        "r7 with 6 N should fail"
    );
}

// Length filter semantics: len < required_length → fail; max_length=0 means no upper bound.
// Default length_required=15.
#[test]
fn se_length_filter_pass_fail() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("out.fq");
    let input = fixture("se_length.fastq");

    run_filter(&[
        "-i",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "-Q", // disable quality filter; test length only
    ]);

    let ids = read_fastq_ids(&out);
    // r1: 20 bp ≥ 15 → PASS
    assert!(
        ids.contains(&"r1_long_enough".to_owned()),
        "r1 20 bp should pass"
    );
    // r2: 8 bp < 15 → FAIL
    assert!(
        !ids.contains(&"r2_too_short".to_owned()),
        "r2 8 bp should fail"
    );
    // r3: 15 bp exactly; 15 < 15 is false → PASS
    assert!(
        ids.contains(&"r3_exactly_required".to_owned()),
        "r3 at exact minimum should pass"
    );
    // r4: 32 bp, no upper limit → PASS
    assert!(
        ids.contains(&"r4_long_but_no_limit".to_owned()),
        "r4 long read no limit should pass"
    );
}

// length_limit > 0: reads longer than limit are dropped.
#[test]
fn se_length_limit_upper_bound() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("out.fq");
    let input = fixture("se_length.fastq");

    // length_limit=20: r4 (32 bp) should be dropped; r1 (20 bp) passes (20 > 20 is false).
    run_filter(&[
        "-i",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "-Q",
        "--length_limit",
        "20",
        "-l",
        "1",
    ]);

    let ids = read_fastq_ids(&out);
    assert!(
        ids.contains(&"r1_long_enough".to_owned()),
        "r1 20 bp at limit should pass"
    );
    assert!(
        !ids.contains(&"r4_long_but_no_limit".to_owned()),
        "r4 32 bp over limit should fail"
    );
}

// Lowercase 'n' must count toward n_base_limit, same as uppercase 'N'.
// Fixture: 20 bp reads with 6 n/N bases; n_base_limit=5 → fail (6 > 5).
// r4 has 1 lowercase n → 1 > 5 is false → pass.
#[test]
fn se_lowercase_n_counts_toward_limit() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("out.fq");
    let input = fixture("se_n_lower.fastq");

    run_filter(&[
        "-i",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "-L",
        "--n_base_limit",
        "5",
    ]);

    let ids = read_fastq_ids(&out);
    assert!(
        !ids.contains(&"r1_uppercase_n_fail".to_owned()),
        "r1 with 6 uppercase N should fail"
    );
    assert!(
        !ids.contains(&"r2_lowercase_n_fail".to_owned()),
        "r2 with 6 lowercase n should fail"
    );
    assert!(
        !ids.contains(&"r3_mixed_case_n_fail".to_owned()),
        "r3 with 6 mixed N/n should fail"
    );
    assert!(
        ids.contains(&"r4_lowercase_n_ok".to_owned()),
        "r4 with 1 lowercase n should pass"
    );
}

// PE: pair is emitted only if BOTH mates pass; if either fails the whole pair is dropped.
#[test]
fn pe_filter_both_must_pass() {
    let tmp = tempfile::tempdir().unwrap();
    let out1 = tmp.path().join("out_r1.fq");
    let out2 = tmp.path().join("out_r2.fq");
    let in1 = fixture("pe_filter.fastq.r1");
    let in2 = fixture("pe_filter.fastq.r2");

    // R1/R2 reads are 20 bp; low-qual limit = 8.0; 9 '!' → fail.
    run_filter(&[
        "-i",
        in1.to_str().unwrap(),
        "-I",
        in2.to_str().unwrap(),
        "-o",
        out1.to_str().unwrap(),
        "-O",
        out2.to_str().unwrap(),
        "-L",
    ]);

    let ids1 = read_fastq_ids(&out1);
    let ids2 = read_fastq_ids(&out2);
    // pair1: both pass → in both outputs
    assert!(
        ids1.contains(&"pair1_both_pass".to_owned()),
        "pair1 R1 should be kept"
    );
    assert!(
        ids2.contains(&"pair1_both_pass".to_owned()),
        "pair1 R2 should be kept"
    );
    // pair2: R1 fails (9 low-qual) → whole pair dropped
    assert!(
        !ids1.contains(&"pair2_r1_fails".to_owned()),
        "pair2 R1 dropped (R1 fails)"
    );
    assert!(
        !ids2.contains(&"pair2_r1_fails".to_owned()),
        "pair2 R2 dropped (R1 fails)"
    );
    // pair3: R2 fails → whole pair dropped
    assert!(
        !ids1.contains(&"pair3_r2_fails".to_owned()),
        "pair3 R1 dropped (R2 fails)"
    );
    assert!(
        !ids2.contains(&"pair3_r2_fails".to_owned()),
        "pair3 R2 dropped (R2 fails)"
    );
    // pair4: both fail → dropped
    assert!(
        !ids1.contains(&"pair4_both_fail".to_owned()),
        "pair4 R1 dropped"
    );
    assert!(
        !ids2.contains(&"pair4_both_fail".to_owned()),
        "pair4 R2 dropped"
    );
    // symmetry: R1 and R2 outputs always have the same set of IDs
    assert_eq!(ids1, ids2, "R1 and R2 outputs must have identical read IDs");
}
