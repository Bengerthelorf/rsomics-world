use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use flate2::Compression;
use flate2::write::GzEncoder;
use needletail::parse_fastx_file;
use rsomics_common::{Context, Result, RsomicsError};

use crate::filter::{FilterConfig, FilterResult, classify};
use crate::polyg::{PolyGConfig, find_polyg_3p};
use crate::report::{FastpJsonReport, FilteringResult};
use crate::stats::ReadStats;
use crate::trim::{AdapterConfig, find_adapter_3p};
use crate::umi::{UmiConfig, UmiLoc, extract as umi_extract};

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
            .rs_with_context(|| format!("creating output FASTQ {}", path.display()))?;
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
            Self::Plain(mut w) => w.flush().rs_context("flushing plain output writer")?,
            Self::Gzip(w) => {
                w.finish().rs_context("finishing gzip output stream")?;
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

fn parse_err(prefix: &str, e: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!("{prefix}: {e}"))
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
        .map_err(|e| parse_err(&format!("opening input FASTQ {}", input.display()), e))?;
    let mut writer = FastqWriter::create(output)?;
    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| parse_err("malformed FASTQ record", e))?;
        rec.write(&mut writer, None)
            .map_err(|e| parse_err("writing record to output", e))?;
    }
    writer.finalize()
}

/// Write one FASTQ record from individual id / seq / qual slices. Needed when
/// trimming has produced new (shorter) seq / qual that diverge from needletail's
/// stored record bytes, so we can't use `record.write`.
fn write_record<W: Write>(
    writer: &mut W,
    id: &[u8],
    seq: &[u8],
    qual: &[u8],
) -> std::io::Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(id)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(qual)?;
    writer.write_all(b"\n")
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

/// Stream a single-end FASTQ through optional adapter trimming, then quality
/// / length / N filters, accumulating per-read statistics and writing only
/// passing reads to `output`. Optionally emit a fastp-compatible JSON report
/// to `json_path`.
///
/// # Errors
///
/// Returns `Err` if input parsing, output writing, or JSON serialization fails.
#[allow(clippy::too_many_arguments)]
pub fn process_se(
    input: &Path,
    output: &Path,
    json_path: Option<&Path>,
    cfg: FilterConfig,
    adapter: Option<&AdapterConfig>,
    polyg: Option<PolyGConfig>,
    umi: Option<UmiConfig>,
) -> Result<SeOutcome> {
    if let Some(u) = umi
        && u.loc != UmiLoc::Read1
    {
        return Err(RsomicsError::ConfigError(
            "single-end UMI extraction only supports umi_loc=read1".into(),
        ));
    }
    let mut reader = parse_fastx_file(input)
        .map_err(|e| parse_err(&format!("opening input FASTQ {}", input.display()), e))?;
    let mut writer = FastqWriter::create(output)?;

    let mut pre = ReadStats::default();
    let mut post = ReadStats::default();
    let mut filtering = FilteringResult::default();

    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| parse_err("malformed FASTQ record", e))?;
        let seq = rec.seq();
        let qual = rec
            .qual()
            .ok_or_else(|| RsomicsError::InvalidInput("FASTQ record missing quality".into()))?;
        pre.observe(&seq, qual);

        let (id_buf, off) = if let Some(u) = umi {
            let Some(pair) = umi_extract(rec.id(), &seq, u) else {
                filtering.record(FilterResult::TooShort);
                continue;
            };
            pair
        } else {
            (rec.id().to_vec(), 0)
        };
        let seq_u = &seq[off..];
        let qual_u = &qual[off..];

        let after_polyg = polyg
            .and_then(|pg| find_polyg_3p(seq_u, pg))
            .unwrap_or(seq_u.len());
        let after_adapter = adapter
            .and_then(|ad| find_adapter_3p(&seq_u[..after_polyg], ad))
            .unwrap_or(after_polyg);
        let trim_at = after_adapter;
        let seq_t = &seq_u[..trim_at];
        let qual_t = &qual_u[..trim_at];

        let outcome = classify(seq_t, qual_t, cfg);
        filtering.record(outcome);
        if matches!(outcome, FilterResult::Pass) {
            post.observe(seq_t, qual_t);
            write_record(&mut writer, &id_buf, seq_t, qual_t)
                .rs_context("writing record to output")?;
        }
    }
    writer.finalize()?;

    if let Some(path) = json_path {
        let report = FastpJsonReport::from_se(
            &pre,
            &post,
            FilteringResult {
                passed_filter_reads: filtering.passed_filter_reads,
                low_quality_reads: filtering.low_quality_reads,
                too_many_n_reads: filtering.too_many_n_reads,
                too_short_reads: filtering.too_short_reads,
            },
        );
        write_json_report(&report, path)?;
    }

    Ok(SeOutcome {
        pre_filter: pre,
        post_filter: post,
        filtering,
    })
}

fn write_json_report(report: &FastpJsonReport, path: &Path) -> Result<()> {
    let mut json_writer = BufWriter::new(
        File::create(path)
            .rs_with_context(|| format!("creating JSON report {}", path.display()))?,
    );
    serde_json::to_writer_pretty(&mut json_writer, report)
        .map_err(|e| parse_err("serializing JSON report", e))?;
    json_writer.flush().rs_context("flushing JSON writer")?;
    Ok(())
}

/// Stream a paired-end FASTQ through the same filter / stats / report pipeline
/// as [`process_se`]. A pair is dropped iff either mate fails any filter;
/// pre/post stats are tracked separately for R1 and R2.
///
/// # Errors
///
/// Returns `Err` if input parsing, output writing, JSON serialization fails,
/// or the two input files have a different number of records.
#[allow(clippy::too_many_lines, clippy::too_many_arguments)]
pub fn process_pe(
    in1: &Path,
    in2: &Path,
    out1: &Path,
    out2: &Path,
    json_path: Option<&Path>,
    cfg: FilterConfig,
    adapter: Option<&AdapterConfig>,
    polyg: Option<PolyGConfig>,
    umi: Option<UmiConfig>,
) -> Result<PeOutcome> {
    let mut r1_reader = parse_fastx_file(in1)
        .map_err(|e| parse_err(&format!("opening input R1 {}", in1.display()), e))?;
    let mut r2_reader = parse_fastx_file(in2)
        .map_err(|e| parse_err(&format!("opening input R2 {}", in2.display()), e))?;
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
                let rec1 = rec1.map_err(|e| parse_err("malformed R1 record", e))?;
                let rec2 = rec2.map_err(|e| parse_err("malformed R2 record", e))?;
                let seq1 = rec1.seq();
                let q1 = rec1
                    .qual()
                    .ok_or_else(|| RsomicsError::InvalidInput("R1 missing quality".into()))?;
                let seq2 = rec2.seq();
                let q2 = rec2
                    .qual()
                    .ok_or_else(|| RsomicsError::InvalidInput("R2 missing quality".into()))?;
                pre_r1.observe(&seq1, q1);
                pre_r2.observe(&seq2, q2);

                // Extract UMI from the configured mate. We stamp the umi
                // string onto BOTH mates' read ids so the pair stays
                // identifiable downstream.
                let (id1_buf, id2_buf, off1, off2) = if let Some(u) = umi {
                    let src_id = rec1.id();
                    let donor = if u.loc == UmiLoc::Read1 { &seq1 } else { &seq2 };
                    let Some((new_id, off)) = umi_extract(src_id, donor, u) else {
                        filtering.record(FilterResult::TooShort);
                        continue;
                    };
                    let id2_buf = stamp_umi(rec2.id(), &new_id[src_id.len()..]);
                    let (o1, o2) = if u.loc == UmiLoc::Read1 {
                        (off, 0)
                    } else {
                        (0, off)
                    };
                    (new_id, id2_buf, o1, o2)
                } else {
                    (rec1.id().to_vec(), rec2.id().to_vec(), 0, 0)
                };
                let seq1_u = &seq1[off1..];
                let q1_u = &q1[off1..];
                let seq2_u = &seq2[off2..];
                let q2_u = &q2[off2..];

                let g1 = polyg
                    .and_then(|pg| find_polyg_3p(seq1_u, pg))
                    .unwrap_or(seq1_u.len());
                let g2 = polyg
                    .and_then(|pg| find_polyg_3p(seq2_u, pg))
                    .unwrap_or(seq2_u.len());
                let t1 = adapter
                    .and_then(|ad| find_adapter_3p(&seq1_u[..g1], ad))
                    .unwrap_or(g1);
                let t2 = adapter
                    .and_then(|ad| find_adapter_3p(&seq2_u[..g2], ad))
                    .unwrap_or(g2);
                let seq1_t = &seq1_u[..t1];
                let q1_t = &q1_u[..t1];
                let seq2_t = &seq2_u[..t2];
                let q2_t = &q2_u[..t2];

                let v1 = classify(seq1_t, q1_t, cfg);
                let v2 = classify(seq2_t, q2_t, cfg);
                let pair_verdict = pair_filter_result(v1, v2);
                filtering.record(pair_verdict);
                if matches!(pair_verdict, FilterResult::Pass) {
                    post_r1.observe(seq1_t, q1_t);
                    post_r2.observe(seq2_t, q2_t);
                    write_record(&mut w1, &id1_buf, seq1_t, q1_t)
                        .rs_context("writing R1 record")?;
                    write_record(&mut w2, &id2_buf, seq2_t, q2_t)
                        .rs_context("writing R2 record")?;
                }
            }
            (None, None) => break,
            (Some(_), None) | (None, Some(_)) => {
                return Err(RsomicsError::InvalidInput(format!(
                    "paired-end inputs have different record counts: {} vs {}",
                    in1.display(),
                    in2.display(),
                )));
            }
        }
    }

    w1.finalize()?;
    w2.finalize()?;

    if let Some(path) = json_path {
        let report = FastpJsonReport::from_pe(
            &pre_r1,
            &post_r1,
            &pre_r2,
            &post_r2,
            FilteringResult {
                passed_filter_reads: filtering.passed_filter_reads,
                low_quality_reads: filtering.low_quality_reads,
                too_many_n_reads: filtering.too_many_n_reads,
                too_short_reads: filtering.too_short_reads,
            },
        );
        write_json_report(&report, path)?;
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

/// Append the already-formed `:UMI` suffix from the donor mate's new id
/// onto a second mate's id, so both mates carry the same umi tag.
fn stamp_umi(id: &[u8], suffix: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(id.len() + suffix.len());
    out.extend_from_slice(id);
    out.extend_from_slice(suffix);
    out
}
