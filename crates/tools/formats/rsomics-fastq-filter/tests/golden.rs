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
    assert!(
        ids.contains(&"r1_all_high_qual".to_owned()),
        "r1 should pass"
    );
    assert!(
        !ids.contains(&"r2_too_many_low_qual".to_owned()),
        "r2 should fail"
    );
    assert!(
        ids.contains(&"r3_exactly_at_limit".to_owned()),
        "r3 at boundary should pass"
    );
    assert!(
        !ids.contains(&"r4_one_over_limit".to_owned()),
        "r4 one over limit should fail"
    );
    assert!(
        !ids.contains(&"r5_all_low_qual".to_owned()),
        "r5 should fail"
    );
    assert!(
        ids.contains(&"r6_n_bases_ok".to_owned()),
        "r6 with 1 N should pass"
    );
    assert!(
        !ids.contains(&"r7_too_many_n".to_owned()),
        "r7 with 6 N should fail"
    );
}

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
    assert!(
        ids.contains(&"r1_long_enough".to_owned()),
        "r1 20 bp should pass"
    );
    assert!(
        !ids.contains(&"r2_too_short".to_owned()),
        "r2 8 bp should fail"
    );
    assert!(
        ids.contains(&"r3_exactly_required".to_owned()),
        "r3 at exact minimum should pass"
    );
    assert!(
        ids.contains(&"r4_long_but_no_limit".to_owned()),
        "r4 long read no limit should pass"
    );
}

#[test]
fn se_length_limit_upper_bound() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("out.fq");
    let input = fixture("se_length.fastq");

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

#[test]
fn pe_filter_both_must_pass() {
    let tmp = tempfile::tempdir().unwrap();
    let out1 = tmp.path().join("out_r1.fq");
    let out2 = tmp.path().join("out_r2.fq");
    let in1 = fixture("pe_filter.fastq.r1");
    let in2 = fixture("pe_filter.fastq.r2");

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
    assert!(
        ids1.contains(&"pair1_both_pass".to_owned()),
        "pair1 R1 should be kept"
    );
    assert!(
        ids2.contains(&"pair1_both_pass".to_owned()),
        "pair1 R2 should be kept"
    );
    assert!(
        !ids1.contains(&"pair2_r1_fails".to_owned()),
        "pair2 R1 dropped (R1 fails)"
    );
    assert!(
        !ids2.contains(&"pair2_r1_fails".to_owned()),
        "pair2 R2 dropped (R1 fails)"
    );
    assert!(
        !ids1.contains(&"pair3_r2_fails".to_owned()),
        "pair3 R1 dropped (R2 fails)"
    );
    assert!(
        !ids2.contains(&"pair3_r2_fails".to_owned()),
        "pair3 R2 dropped (R2 fails)"
    );
    assert!(
        !ids1.contains(&"pair4_both_fail".to_owned()),
        "pair4 R1 dropped"
    );
    assert!(
        !ids2.contains(&"pair4_both_fail".to_owned()),
        "pair4 R2 dropped"
    );
    assert_eq!(ids1, ids2, "R1 and R2 outputs must have identical read IDs");
}
