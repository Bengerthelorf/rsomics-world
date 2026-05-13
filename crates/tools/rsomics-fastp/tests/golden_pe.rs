//! Paired-end golden-fixture tests.

use std::fs;
use std::path::PathBuf;

use rsomics_fastp::filter::FilterConfig;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("golden")
        .join(name)
}

#[test]
fn process_pe_rejects_pair_when_either_mate_fails() {
    let in1 = fixture("pe_mixed_r1.fastq");
    let in2 = fixture("pe_mixed_r2.fastq");
    let out1 = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .expect("tempfile");
    let out2 = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .expect("tempfile");
    let json = tempfile::Builder::new()
        .suffix(".json")
        .tempfile()
        .expect("tempfile");

    let outcome = rsomics_fastp::io::process_pe(
        &in1,
        &in2,
        out1.path(),
        out2.path(),
        Some(json.path()),
        FilterConfig::default(),
        None,
        None,
    )
    .expect("process_pe");

    // pair_001: both pass. pair_002: R1 too short. pair_003: R2 low quality.
    // pair_004: both pass. → 2 pass, 1 too short, 1 low quality.
    assert_eq!(outcome.filtering.passed_filter_reads, 2);
    assert_eq!(outcome.filtering.too_short_reads, 1);
    assert_eq!(outcome.filtering.low_quality_reads, 1);
    assert_eq!(outcome.filtering.too_many_n_reads, 0);
    assert_eq!(outcome.pre_filter_r1.total_reads, 4);
    assert_eq!(outcome.pre_filter_r2.total_reads, 4);
    assert_eq!(outcome.post_filter_r1.total_reads, 2);
    assert_eq!(outcome.post_filter_r2.total_reads, 2);

    // Output FASTQs contain exactly the surviving pairs in order.
    let r1_text = fs::read_to_string(out1.path()).expect("read r1");
    let r2_text = fs::read_to_string(out2.path()).expect("read r2");
    assert!(r1_text.contains("@pair_001"));
    assert!(r1_text.contains("@pair_004"));
    assert!(!r1_text.contains("@pair_002"));
    assert!(!r1_text.contains("@pair_003"));
    assert!(r2_text.contains("@pair_001"));
    assert!(r2_text.contains("@pair_004"));
    assert!(!r2_text.contains("@pair_002"));
    assert!(!r2_text.contains("@pair_003"));

    // JSON report's aggregate summary counts BOTH mates.
    let json_text = fs::read_to_string(json.path()).expect("read json");
    let parsed: serde_json::Value = serde_json::from_str(&json_text).expect("json parse");
    assert_eq!(parsed["summary"]["before_filtering"]["total_reads"], 8);
    assert_eq!(parsed["summary"]["after_filtering"]["total_reads"], 4);
    assert_eq!(parsed["filtering_result"]["passed_filter_reads"], 2);
}
