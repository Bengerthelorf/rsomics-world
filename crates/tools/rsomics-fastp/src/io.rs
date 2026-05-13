use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use needletail::parse_fastx_file;

use crate::filter::{FilterConfig, FilterResult, classify};
use crate::report::{FastpJsonReport, FilteringResult};
use crate::stats::ReadStats;

/// Identity-copy a single-end FASTQ file. No transformation; validates the
/// reader / writer plumbing in isolation before filtering layers ride on top.
///
/// # Errors
///
/// Returns `Err` if the input cannot be opened, a record fails to parse,
/// the output cannot be created, or a write to the output fails.
pub fn copy_se(input: &Path, output: &Path) -> Result<()> {
    let mut reader = parse_fastx_file(input)
        .with_context(|| format!("opening input FASTQ {}", input.display()))?;
    let mut writer = BufWriter::new(
        File::create(output)
            .with_context(|| format!("creating output FASTQ {}", output.display()))?,
    );
    while let Some(record) = reader.next() {
        let rec = record.context("malformed FASTQ record")?;
        rec.write(&mut writer, None)
            .context("writing record to output")?;
    }
    writer.flush().context("flushing output writer")?;
    Ok(())
}

/// Outcome of a single-end preprocessing run — both pre- and post-filter
/// statistics, plus the per-category filter counts.
#[derive(Debug)]
pub struct SeOutcome {
    pub pre_filter: ReadStats,
    pub post_filter: ReadStats,
    pub filtering: FilteringResult,
}

/// Stream a single-end FASTQ through quality / length / N filters and accumulate
/// per-read statistics, writing only the passing reads to `output`. Optionally
/// emit a fastp-compatible JSON report to `json_path`.
///
/// # Errors
///
/// Returns `Err` if input parsing, output writing, or JSON serialization fails.
pub fn process_se(
    input: &Path,
    output: &Path,
    json_path: Option<&Path>,
    cfg: FilterConfig,
) -> Result<SeOutcome> {
    let mut reader = parse_fastx_file(input)
        .with_context(|| format!("opening input FASTQ {}", input.display()))?;
    let mut writer = BufWriter::new(
        File::create(output)
            .with_context(|| format!("creating output FASTQ {}", output.display()))?,
    );

    let mut pre = ReadStats::default();
    let mut post = ReadStats::default();
    let mut filtering = FilteringResult::default();

    while let Some(record) = reader.next() {
        let rec = record.context("malformed FASTQ record")?;
        let seq = rec.seq();
        let qual = rec.qual().context("FASTQ record missing quality scores")?;
        pre.observe(&seq, qual);

        let outcome = classify(&seq, qual, cfg);
        filtering.record(outcome);
        if matches!(outcome, FilterResult::Pass) {
            post.observe(&seq, qual);
            rec.write(&mut writer, None)
                .context("writing record to output")?;
        }
    }
    writer.flush().context("flushing output writer")?;

    if let Some(path) = json_path {
        let report = FastpJsonReport::from_stats(
            &pre,
            &post,
            FilteringResult {
                passed_filter_reads: filtering.passed_filter_reads,
                low_quality_reads: filtering.low_quality_reads,
                too_many_n_reads: filtering.too_many_n_reads,
                too_short_reads: filtering.too_short_reads,
            },
        );
        let mut json_writer = BufWriter::new(
            File::create(path)
                .with_context(|| format!("creating JSON report {}", path.display()))?,
        );
        serde_json::to_writer_pretty(&mut json_writer, &report)
            .context("serializing JSON report")?;
        json_writer.flush().context("flushing JSON writer")?;
    }

    Ok(SeOutcome {
        pre_filter: pre,
        post_filter: post,
        filtering,
    })
}
