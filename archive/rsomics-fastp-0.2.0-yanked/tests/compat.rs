//! Behavioural-compatibility tests against upstream `fastp`.
//!
//! Two tiers:
//!
//! - **Tier 1** (always-on): small synthetic SE FASTQ committed under
//!   `tests/golden/`. Asserts surviving-read and surviving-base counts
//!   match upstream when `fastp` is on PATH.
//! - **Tier 2** (download-gated): chr22 wgsim PE fixture from
//!   `tests/fixtures-manifest.toml`. Asserts JSON summary scalars
//!   (`total_reads`, `passed_filter_reads`) agree, plus output FASTQ
//!   surviving record / base counts agree.
//!
//! CI installs `fastp` so both tiers run there; locally they skip
//! gracefully if the upstream binary or network is unavailable.

use std::path::Path;
use std::process::Command;

use rsomics_common::test_support::{tier2, tool_on_path};
use rsomics_common::{fixture_path, tier2_manifest_path};

fn count_fastq_records(path: &Path) -> usize {
    let text = std::fs::read_to_string(path).expect("read fastq");
    text.lines().filter(|l| l.starts_with('@')).count()
}

fn sum_fastq_bases(path: &Path) -> usize {
    let text = std::fs::read_to_string(path).expect("read fastq");
    // Every 4-line FASTQ block is: @id / seq / + / qual. The second line of
    // each block is the sequence we count.
    let lines: Vec<&str> = text.lines().collect();
    lines
        .chunks_exact(4)
        .map(|chunk| chunk[1].len())
        .sum::<usize>()
}

fn run_upstream_fastp(input: &Path, output: &Path) {
    let status = Command::new("fastp")
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(output)
        .arg("--disable_adapter_trimming")
        .arg("--json")
        .arg(output.with_extension("fastp.json"))
        .arg("--html")
        .arg(output.with_extension("fastp.html"))
        .status()
        .expect("spawn fastp");
    assert!(status.success(), "upstream fastp exited non-zero");
}

#[test]
fn surviving_read_count_matches_upstream_on_se_mixed() {
    if !tool_on_path("fastp") {
        eprintln!("skip: upstream fastp not on PATH");
        return;
    }
    let input = fixture_path!("se_mixed.fastq");
    let ours = tempfile::Builder::new()
        .suffix("_ours.fastq")
        .tempfile()
        .expect("tempfile");
    let theirs = tempfile::Builder::new()
        .suffix("_theirs.fastq")
        .tempfile()
        .expect("tempfile");

    rsomics_fastp::io::process_se(
        &input,
        ours.path(),
        None,
        rsomics_fastp::filter::FilterConfig::default(),
        None,
        None,
        None,
    )
    .expect("rsomics-fastp run");
    run_upstream_fastp(&input, theirs.path());

    let ours_n = count_fastq_records(ours.path());
    let theirs_n = count_fastq_records(theirs.path());
    assert_eq!(
        ours_n, theirs_n,
        "surviving read count differs: ours={ours_n}, fastp={theirs_n}"
    );

    let ours_b = sum_fastq_bases(ours.path());
    let theirs_b = sum_fastq_bases(theirs.path());
    assert_eq!(
        ours_b, theirs_b,
        "surviving base count differs: ours={ours_b}, fastp={theirs_b}"
    );
}

fn count_records_gz(path: &Path) -> usize {
    use std::io::Read;
    let f = std::fs::File::open(path).expect("open gz");
    let mut decoder = flate2::read::MultiGzDecoder::new(f);
    let mut buf = String::new();
    decoder.read_to_string(&mut buf).expect("decode gz");
    buf.lines().filter(|l| l.starts_with('@')).count()
}

fn sum_bases_gz(path: &Path) -> usize {
    use std::io::Read;
    let f = std::fs::File::open(path).expect("open gz");
    let mut decoder = flate2::read::MultiGzDecoder::new(f);
    let mut buf = String::new();
    decoder.read_to_string(&mut buf).expect("decode gz");
    buf.lines()
        .collect::<Vec<_>>()
        .chunks_exact(4)
        .map(|c| c[1].len())
        .sum()
}

#[test]
fn pe_surviving_pair_and_base_counts_match_upstream_on_chr22_tier2() {
    if !tool_on_path("fastp") {
        eprintln!("skip: upstream fastp not on PATH");
        return;
    }
    let manifest = tier2_manifest_path!();
    let r1 = match tier2::fetch(&manifest, "chr22_sub_r1.fastq.gz") {
        Ok(p) => p,
        Err(e) => {
            eprintln!("skip: tier-2 R1 fetch failed ({e})");
            return;
        }
    };
    let r2 = match tier2::fetch(&manifest, "chr22_sub_r2.fastq.gz") {
        Ok(p) => p,
        Err(e) => {
            eprintln!("skip: tier-2 R2 fetch failed ({e})");
            return;
        }
    };

    let tmpdir = tempfile::tempdir().expect("tempdir");
    let our_r1 = tmpdir.path().join("our_r1.fastq.gz");
    let our_r2 = tmpdir.path().join("our_r2.fastq.gz");
    let our_json = tmpdir.path().join("our.json");
    let their_r1 = tmpdir.path().join("their_r1.fastq.gz");
    let their_r2 = tmpdir.path().join("their_r2.fastq.gz");
    let their_json = tmpdir.path().join("their.json");

    rsomics_fastp::io::process_pe(
        &r1,
        &r2,
        &our_r1,
        &our_r2,
        Some(&our_json),
        rsomics_fastp::filter::FilterConfig::default(),
        None,
        None,
        None,
    )
    .expect("rsomics-fastp PE run");

    let status = Command::new("fastp")
        .arg("-i")
        .arg(&r1)
        .arg("-I")
        .arg(&r2)
        .arg("-o")
        .arg(&their_r1)
        .arg("-O")
        .arg(&their_r2)
        .arg("--disable_adapter_trimming")
        .arg("--json")
        .arg(&their_json)
        .arg("--html")
        .arg(tmpdir.path().join("their.html"))
        .arg("--thread")
        .arg("1")
        .stderr(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .expect("spawn fastp");
    assert!(status.success(), "upstream fastp exited non-zero");

    let our_n = count_records_gz(&our_r1) + count_records_gz(&our_r2);
    let their_n = count_records_gz(&their_r1) + count_records_gz(&their_r2);
    assert_eq!(
        our_n, their_n,
        "surviving record count (R1+R2) differs: ours={our_n}, fastp={their_n}"
    );

    let our_b = sum_bases_gz(&our_r1) + sum_bases_gz(&our_r2);
    let their_b = sum_bases_gz(&their_r1) + sum_bases_gz(&their_r2);
    assert_eq!(
        our_b, their_b,
        "surviving base count (R1+R2) differs: ours={our_b}, fastp={their_b}"
    );

    // Cross-check JSON scalars. Both reports must agree on the headline
    // numbers; per-cycle curves are allowed to diverge (fastp's are real,
    // ours mirror the mean — documented).
    let rsomics: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&our_json).expect("read our json"))
            .expect("parse our json");
    let upstream: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&their_json).expect("read their json"))
            .expect("parse their json");
    assert_eq!(
        rsomics["summary"]["after_filtering"]["total_reads"],
        upstream["summary"]["after_filtering"]["total_reads"],
        "JSON summary.after_filtering.total_reads mismatch"
    );
    assert_eq!(
        rsomics["filtering_result"]["passed_filter_reads"],
        upstream["filtering_result"]["passed_filter_reads"],
        "JSON filtering_result.passed_filter_reads mismatch"
    );
}
