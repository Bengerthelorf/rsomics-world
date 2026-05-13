//! Golden-fixture tests that don't depend on the upstream fastp binary.
//! These exercise the internal pipeline end-to-end against small hand-curated
//! FASTQ files under `tests/golden/`.

use std::fs;
use std::path::PathBuf;

use rsomics_fastp::filter::FilterConfig;
use rsomics_fastp::trim::AdapterConfig;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("golden")
        .join(name)
}

#[test]
fn copy_se_is_identity() {
    let input = fixture("se_basic.fastq");
    let tmp = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .expect("tempfile");

    rsomics_fastp::io::copy_se(&input, tmp.path()).expect("copy_se");

    let expected = fs::read(&input).expect("read input");
    let actual = fs::read(tmp.path()).expect("read output");
    assert_eq!(
        actual, expected,
        "identity copy must preserve the input byte-for-byte"
    );
}

#[test]
fn copy_se_writes_gzipped_output_when_extension_is_gz() {
    use std::io::Read;

    let input = fixture("se_basic.fastq");
    let tmp = tempfile::Builder::new()
        .suffix(".fastq.gz")
        .tempfile()
        .expect("tempfile");

    rsomics_fastp::io::copy_se(&input, tmp.path()).expect("copy_se");

    // Output must start with the gzip magic bytes.
    let mut head = [0u8; 2];
    std::fs::File::open(tmp.path())
        .expect("open gz output")
        .read_exact(&mut head)
        .expect("read magic");
    assert_eq!(head, [0x1f, 0x8b], "gzipped output missing gzip magic");

    // Round-tripping through needletail must recover the original record content.
    let mut reader = needletail::parse_fastx_file(tmp.path()).expect("parse gz");
    let expected = std::fs::read(&input).expect("read input");
    let mut decoded = Vec::with_capacity(expected.len());
    while let Some(rec) = reader.next() {
        rec.expect("record")
            .write(&mut decoded, None)
            .expect("write");
    }
    assert_eq!(
        decoded, expected,
        "gzipped round-trip must recover the original byte stream"
    );
}

#[test]
fn process_se_classifies_each_failure_mode() {
    let input = fixture("se_mixed.fastq");
    let out = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .expect("tempfile");
    let json = tempfile::Builder::new()
        .suffix(".json")
        .tempfile()
        .expect("tempfile");

    let outcome = rsomics_fastp::io::process_se(
        &input,
        out.path(),
        Some(json.path()),
        FilterConfig::default(),
        None,
    )
    .expect("process_se");

    // Fixture covers each failure mode once + one pass.
    assert_eq!(outcome.filtering.passed_filter_reads, 1);
    assert_eq!(outcome.filtering.too_short_reads, 1);
    assert_eq!(outcome.filtering.too_many_n_reads, 1);
    assert_eq!(outcome.filtering.low_quality_reads, 1);
    assert_eq!(outcome.pre_filter.total_reads, 4);
    assert_eq!(outcome.post_filter.total_reads, 1);

    // JSON report renders and parses back.
    let json_text = fs::read_to_string(json.path()).expect("read json");
    let parsed: serde_json::Value = serde_json::from_str(&json_text).expect("json parse");
    assert_eq!(parsed["filtering_result"]["passed_filter_reads"], 1);
    assert_eq!(parsed["filtering_result"]["too_many_N_reads"], 1);
    assert_eq!(parsed["summary"]["before_filtering"]["total_reads"], 4);
    assert_eq!(parsed["summary"]["after_filtering"]["total_reads"], 1);

    // Only the passing read should appear in the output FASTQ.
    let out_text = fs::read_to_string(out.path()).expect("read output");
    assert!(out_text.contains("@pass_high_q"));
    assert!(!out_text.contains("@too_short"));
    assert!(!out_text.contains("@too_many_n"));
    assert!(!out_text.contains("@low_quality"));
}

#[test]
fn process_se_trims_adapter_at_3prime() {
    let input = fixture("se_adapter.fastq");
    let out = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .expect("tempfile");

    let outcome = rsomics_fastp::io::process_se(
        &input,
        out.path(),
        None,
        FilterConfig::default(),
        Some(&AdapterConfig::illumina_truseq_r1()),
    )
    .expect("process_se");

    // Fixture has 2 reads, both with the TruSeq adapter at offset 20.
    // After trimming both should pass (length ≥ 15) and the post stats should
    // show 20 bp insert each (40 total bases).
    assert_eq!(outcome.filtering.passed_filter_reads, 2);
    assert_eq!(outcome.post_filter.total_reads, 2);
    assert_eq!(outcome.post_filter.total_bases, 40);

    // Output FASTQ records should be the 20-bp inserts, no adapter bytes.
    let out_text = fs::read_to_string(out.path()).expect("read output");
    assert!(out_text.contains("ACGTACGTACGTACGTACGT\n"));
    assert!(
        !out_text.contains("AGATCGGAAGAG"),
        "adapter bytes leaked into output"
    );
}
