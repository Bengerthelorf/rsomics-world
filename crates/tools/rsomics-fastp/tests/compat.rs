//! Behavioural-compatibility tests against upstream `fastp`.
//!
//! These run `rsomics-fastp` and upstream `fastp` on the same golden FASTQ
//! and compare invariants we expect to hold across implementations
//! (filtered-read count, surviving-base count). When the upstream binary
//! isn't on PATH the tests print a skip notice and pass — CI is responsible
//! for ensuring fastp is installed where compat is meaningful.

use std::path::Path;
use std::process::Command;

use rsomics_common::{fixture_path, test_support::tool_on_path};

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
