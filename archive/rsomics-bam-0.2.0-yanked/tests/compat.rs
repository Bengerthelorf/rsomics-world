//! Behavioural-compatibility tests for `rsomics-bam` against upstream
//! `samtools`.
//!
//! Two tiers of coverage:
//!
//! - **Tier 1** (always-on): synthetic 17-record BAM built at test setup
//!   via [`mod synth`]. Verifies the library API (`count_records`) and
//!   compares against `samtools view -c` when samtools is on PATH.
//! - **Tier 2** (download-gated): real-shape 100k-record BAM from chr22
//!   simulated reads, fetched per `tests/fixtures-manifest.toml`. Skipped
//!   when network or samtools is unavailable.
//!
//! Lives in `tests/compat.rs` per CLAUDE.md's Layer-B requirement.

mod synth;

use std::process::Command;

use rsomics_common::test_support::{tier2, tool_on_path};
use rsomics_common::tier2_manifest_path;

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

#[test]
fn rsomics_bam_view_count_matches_samtools_on_chr22_tier2() {
    if !tool_on_path("samtools") {
        eprintln!("skip: samtools not on PATH");
        return;
    }
    let manifest = tier2_manifest_path!();
    let bam = match tier2::fetch(&manifest, "chr22_sub.bam") {
        Ok(p) => p,
        Err(e) => {
            eprintln!("skip: tier-2 fixture fetch failed ({e})");
            return;
        }
    };

    let ours = rsomics_bam::cmd::view::count_records(&bam).expect("rsomics-bam count_records");

    let out = Command::new("samtools")
        .arg("view")
        .arg("-c")
        .arg(&bam)
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
        "rsomics-bam count {ours} != samtools view -c {theirs} on chr22 tier-2 fixture"
    );
    // The fixture is the 50k-paired-end wgsim simulation; minimap2-sr
    // emits one primary alignment per mate ≈ 100k records.
    assert_eq!(
        ours, 100_000,
        "fixture expected 100k records, got {ours} — manifest may be stale"
    );
}
