use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;
use needletail::parse_fastx_file;

use crate::filter::{FilterConfig, FilterResult, classify};
use crate::report::{FastpJsonReport, FilteringResult};
use crate::stats::ReadStats;

/// FASTQ output sink. Both arms write through a `BufWriter` so needletail's
/// small per-record writes batch into larger I/O. The gzip variant must be
/// `finalize`d to emit the gzip trailer cleanly; `Drop` calls `try_finish`
/// which writes the trailer but silently swallows late errors (e.g. disk full
/// during the final flush), so the explicit `finalize` is the supported path.
enum FastqWriter {
    Plain(BufWriter<File>),
    Gzip(GzEncoder<BufWriter<File>>),
}

impl FastqWriter {
    fn create(path: &Path) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("creating output FASTQ {}", path.display()))?;
        let buf = BufWriter::new(file);
        if path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("gz"))
        {
            Ok(Self::Gzip(GzEncoder::new(buf, Compression::default())))
        } else {
            Ok(Self::Plain(buf))
        }
    }

    fn finalize(self) -> Result<()> {
        match self {
            Self::Plain(mut w) => w.flush().context("flushing plain output writer")?,
            Self::Gzip(w) => {
                w.finish().context("finishing gzip output stream")?;
            }
        }
        Ok(())
    }
}

impl Write for FastqWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Plain(w) => w.write(buf),
            Self::Gzip(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Plain(w) => w.flush(),
            Self::Gzip(w) => w.flush(),
        }
    }
}

/// Identity-copy a single-end FASTQ file. No transformation; validates the
/// reader / writer plumbing in isolation before filtering layers ride on top.
/// Input compression is auto-detected by needletail; output is gzipped iff the
/// path ends in `.gz`.
///
/// # Errors
///
/// Returns `Err` if the input cannot be opened, a record fails to parse,
/// the output cannot be created, or a write to the output fails.
pub fn copy_se(input: &Path, output: &Path) -> Result<()> {
    let mut reader = parse_fastx_file(input)
        .with_context(|| format!("opening input FASTQ {}", input.display()))?;
    let mut writer = FastqWriter::create(output)?;
    while let Some(record) = reader.next() {
        let rec = record.context("malformed FASTQ record")?;
        rec.write(&mut writer, None)
            .context("writing record to output")?;
    }
    writer.finalize()
}

/// Outcome of a single-end preprocessing run — both pre- and post-filter
/// statistics, plus the per-category filter counts.
#[derive(Debug)]
pub struct SeOutcome {
    pub pre_filter: ReadStats,
    pub post_filter: ReadStats,
    pub filtering: FilteringResult,
}

/// Outcome of a paired-end preprocessing run. `filtering` counts are at the
/// pair level (a pair is rejected as soon as either mate fails any filter);
/// pre/post stats are tracked separately for R1 and R2 plus an aggregate.
#[derive(Debug)]
pub struct PeOutcome {
    pub pre_filter_r1: ReadStats,
    pub pre_filter_r2: ReadStats,
    pub post_filter_r1: ReadStats,
    pub post_filter_r2: ReadStats,
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
    let mut writer = FastqWriter::create(output)?;

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
    writer.finalize()?;

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

/// Stream a paired-end FASTQ through the same filter / stats / report pipeline
/// as [`process_se`]. A pair is dropped iff either mate fails any filter;
/// pre/post stats are tracked separately for R1 and R2.
///
/// # Errors
///
/// Returns `Err` if input parsing, output writing, JSON serialization fails,
/// or the two input files have a different number of records.
#[allow(clippy::too_many_lines)]
pub fn process_pe(
    in1: &Path,
    in2: &Path,
    out1: &Path,
    out2: &Path,
    json_path: Option<&Path>,
    cfg: FilterConfig,
) -> Result<PeOutcome> {
    let mut r1_reader =
        parse_fastx_file(in1).with_context(|| format!("opening input R1 {}", in1.display()))?;
    let mut r2_reader =
        parse_fastx_file(in2).with_context(|| format!("opening input R2 {}", in2.display()))?;
    let mut w1 = FastqWriter::create(out1)?;
    let mut w2 = FastqWriter::create(out2)?;

    let mut pre_r1 = ReadStats::default();
    let mut pre_r2 = ReadStats::default();
    let mut post_r1 = ReadStats::default();
    let mut post_r2 = ReadStats::default();
    let mut filtering = FilteringResult::default();

    loop {
        let r1 = r1_reader.next();
        let r2 = r2_reader.next();
        match (r1, r2) {
            (Some(rec1), Some(rec2)) => {
                let rec1 = rec1.context("malformed R1 record")?;
                let rec2 = rec2.context("malformed R2 record")?;
                let seq1 = rec1.seq();
                let q1 = rec1.qual().context("R1 missing quality")?;
                let seq2 = rec2.seq();
                let q2 = rec2.qual().context("R2 missing quality")?;
                pre_r1.observe(&seq1, q1);
                pre_r2.observe(&seq2, q2);

                let v1 = classify(&seq1, q1, cfg);
                let v2 = classify(&seq2, q2, cfg);
                let pair_verdict = pair_filter_result(v1, v2);
                filtering.record(pair_verdict);
                if matches!(pair_verdict, FilterResult::Pass) {
                    post_r1.observe(&seq1, q1);
                    post_r2.observe(&seq2, q2);
                    rec1.write(&mut w1, None).context("writing R1 record")?;
                    rec2.write(&mut w2, None).context("writing R2 record")?;
                }
            }
            (None, None) => break,
            (Some(_), None) | (None, Some(_)) => {
                anyhow::bail!(
                    "paired-end inputs have different record counts: {} vs {}",
                    in1.display(),
                    in2.display(),
                );
            }
        }
    }

    w1.finalize()?;
    w2.finalize()?;

    if let Some(path) = json_path {
        let mut pre_agg = pre_r1.clone();
        merge_stats(&mut pre_agg, &pre_r2);
        let mut post_agg = post_r1.clone();
        merge_stats(&mut post_agg, &post_r2);
        let report = FastpJsonReport::from_stats(
            &pre_agg,
            &post_agg,
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

    Ok(PeOutcome {
        pre_filter_r1: pre_r1,
        pre_filter_r2: pre_r2,
        post_filter_r1: post_r1,
        post_filter_r2: post_r2,
        filtering,
    })
}

/// Reduce two per-mate verdicts to a single pair-level verdict. A pair passes
/// iff both mates pass; otherwise the first failure encountered (R1 before R2)
/// is reported, with the precedence inside each verdict preserved by
/// [`classify`].
fn pair_filter_result(v1: FilterResult, v2: FilterResult) -> FilterResult {
    if matches!(v1, FilterResult::Pass) {
        v2
    } else {
        v1
    }
}

/// Fold the contents of `other` into `into`, summing counts. Used to build
/// the aggregate pre/post stats for the JSON `summary` section from per-mate
/// stats.
fn merge_stats(into: &mut ReadStats, other: &ReadStats) {
    into.total_reads += other.total_reads;
    into.total_bases += other.total_bases;
    into.q20_bases += other.q20_bases;
    into.q30_bases += other.q30_bases;
    into.gc_bases += other.gc_bases;
    into.n_bases += other.n_bases;
}
