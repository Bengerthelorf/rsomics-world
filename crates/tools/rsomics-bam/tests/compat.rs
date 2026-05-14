//! Golden-fixture test for `rsomics-bam view -c`. Builds a small BAM at
//! test setup via [`mod synth`], counts records two ways:
//!
//! - directly via `rsomics_bam::cmd::view::count_records` (library call)
//! - if `samtools` is on PATH, also `samtools view -c <file>` as compat
//!
//! Both must agree on the count.

mod synth;

use std::process::Command;

use rsomics_common::test_support::tool_on_path;

const N: usize = 17;

#[test]
fn rsomics_bam_view_count_matches_record_count() {
    let tmp = tempfile::Builder::new()
        .suffix(".bam")
        .tempfile()
        .expect("tempfile");
    let written = synth::write_unaligned_bam(tmp.path(), N);
    assert_eq!(
        written, N,
        "fixture builder should report N records written"
    );

    let counted = rsomics_bam::cmd::view::count_records(tmp.path()).expect("count");
    assert_eq!(usize::try_from(counted).expect("count fits in usize"), N);
}

#[test]
fn rsomics_bam_view_count_matches_samtools_view_c() {
    if !tool_on_path("samtools") {
        eprintln!("skip: samtools not on PATH");
        return;
    }
    let tmp = tempfile::Builder::new()
        .suffix(".bam")
        .tempfile()
        .expect("tempfile");
    let _ = synth::write_unaligned_bam(tmp.path(), N);

    let ours = rsomics_bam::cmd::view::count_records(tmp.path()).expect("count");

    let out = Command::new("samtools")
        .arg("view")
        .arg("-c")
        .arg(tmp.path())
        .output()
        .expect("spawn samtools");
    assert!(out.status.success(), "samtools view -c exited non-zero");
    let theirs: u64 = String::from_utf8(out.stdout)
        .expect("samtools stdout utf8")
        .trim()
        .parse()
        .expect("samtools count parse");

    assert_eq!(
        ours, theirs,
        "rsomics-bam count {ours} != samtools view -c {theirs}"
    );
}
