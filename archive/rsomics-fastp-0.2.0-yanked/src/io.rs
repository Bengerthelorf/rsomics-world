use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use flate2::Compression;
use flate2::write::GzEncoder;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use rsomics_common::{Context, Result, RsomicsError};

use crate::filter::{FilterConfig, FilterResult, classify};
use crate::polyg::{PolyGConfig, find_polyg_3p};
use crate::report::{FastpJsonReport, FilteringResult};
use crate::stats::ReadStats;
use crate::trim::{AdapterConfig, find_adapter_3p};
use crate::umi::{UmiConfig, UmiLoc, extract as umi_extract};

/// Chunk size for the parallel scatter/gather pipeline. Trade-off:
/// larger chunks amortise the rayon dispatch overhead, smaller chunks
/// reduce memory peak (each in-flight record holds id+seq+qual+stats).
/// 1024 records × ~150 bp ≈ 1.5 MB seq bytes per chunk; comfortable for
/// any modern machine.
const CHUNK_RECORDS: usize = 1024;

/// One FASTQ record decoupled from needletail's borrowed-buffer reader.
/// Owned bytes so the chunk can be processed in parallel without holding
/// the reader's lifetime.
struct OwnedSeRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

/// Result of processing one SE record on a worker thread. The local
/// `ReadStats` carry just this record's contribution; the gather phase
/// folds them into the run-wide totals via `ReadStats::merge`. Owned
/// buffers so the worker can return them across the rayon boundary.
struct ProcessedSe {
    pre_local: ReadStats,
    post_local: ReadStats,
    verdict: FilterResult,
    write: Option<(Vec<u8>, Vec<u8>, Vec<u8>)>,
}

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

/// Pure per-record transform: UMI extract → poly-G trim → adapter trim
/// → quality / length / N classify → optional output buffer build, plus
/// per-record `ReadStats` contributions. Free of any I/O so workers in a
/// rayon pool can run it without contention on the reader or writer.
///
/// Returns `None` from the UMI branch if the read is shorter than the
/// requested UMI length — callers treat that as a `TooShort` outcome
/// for accounting.
///
/// Takes `rec` by value so rayon's `par_iter().map(...)` can pass each
/// owned record into the closure without lifetime gymnastics. The
/// fields are consumed via clones / slice copies, so passing by reference
/// would just push the consumption decision to the caller.
#[allow(clippy::needless_pass_by_value)]
fn process_se_record(
    rec: OwnedSeRecord,
    cfg: FilterConfig,
    adapter: Option<&AdapterConfig>,
    polyg: Option<PolyGConfig>,
    umi: Option<UmiConfig>,
) -> ProcessedSe {
    let mut pre_local = ReadStats::default();
    pre_local.observe(&rec.seq, &rec.qual);

    let (id_buf, off) = if let Some(u) = umi {
        match umi_extract(&rec.id, &rec.seq, u) {
            Some(pair) => pair,
            None => {
                return ProcessedSe {
                    pre_local,
                    post_local: ReadStats::default(),
                    verdict: FilterResult::TooShort,
                    write: None,
                };
            }
        }
    } else {
        (rec.id.clone(), 0)
    };
    let seq_u = &rec.seq[off..];
    let qual_u = &rec.qual[off..];

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
    let (post_local, write) = if matches!(outcome, FilterResult::Pass) {
        let mut p = ReadStats::default();
        p.observe(seq_t, qual_t);
        (p, Some((id_buf, seq_t.to_vec(), qual_t.to_vec())))
    } else {
        (ReadStats::default(), None)
    };
    ProcessedSe {
        pre_local,
        post_local,
        verdict: outcome,
        write,
    }
}

/// Stream a single-end FASTQ through optional UMI / poly-G / adapter
/// trimming and a quality / length / N filter. Records are read in
/// chunks of [`CHUNK_RECORDS`], scattered across the global rayon pool
/// for the CPU-heavy per-record processing, then gathered serially for
/// stats merging and ordered output. `--threads N` is honoured via the
/// pool size set by [`rsomics_common::CommonFlags::install_rayon_pool`].
///
/// Optionally emits a fastp-compatible JSON report to `json_path`.
///
/// # Errors
///
/// Returns `Err` if input parsing, output writing, or JSON serialization
/// fails. UMI-too-short reads (when `umi` is set and the read is shorter
/// than the requested UMI length) count as `TooShort` rather than error.
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
    let mut chunk: Vec<OwnedSeRecord> = Vec::with_capacity(CHUNK_RECORDS);

    loop {
        chunk.clear();
        while chunk.len() < CHUNK_RECORDS {
            match reader.next() {
                Some(r) => {
                    let rec = r.map_err(|e| parse_err("malformed FASTQ record", e))?;
                    let qual = rec.qual().ok_or_else(|| {
                        RsomicsError::InvalidInput("FASTQ record missing quality".into())
                    })?;
                    chunk.push(OwnedSeRecord {
                        id: rec.id().to_vec(),
                        seq: rec.seq().into_owned(),
                        qual: qual.to_vec(),
                    });
                }
                None => break,
            }
        }
        if chunk.is_empty() {
            break;
        }

        // Scatter: parallel transform per record. par_drain preserves
        // order in the collected output, so the writer below sees records
        // in the same order the reader produced them.
        let processed: Vec<ProcessedSe> = chunk
            .par_drain(..)
            .map(|rec| process_se_record(rec, cfg, adapter, polyg, umi))
            .collect();

        // Gather: serial merge of per-record stats + ordered writes.
        for p in processed {
            pre.merge(&p.pre_local);
            filtering.record(p.verdict);
            if let Some((id, seq, qual)) = p.write {
                post.merge(&p.post_local);
                write_record(&mut writer, &id, &seq, &qual).rs_context("writing record")?;
            }
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

/// One paired-end record decoupled from needletail's readers. Owned bytes
/// so chunks can scatter across rayon workers without holding either
/// reader's lifetime.
struct OwnedPeRecord {
    rec1: OwnedSeRecord,
    rec2: OwnedSeRecord,
}

/// Triple of (id, seq, qual) for one writable record, post-trim.
type WriteRec = (Vec<u8>, Vec<u8>, Vec<u8>);

/// Per-record worker output for PE. The `write` buffers are pre-sliced
/// to the post-trim region so the writer just emits them verbatim.
struct ProcessedPe {
    pre_local_r1: ReadStats,
    pre_local_r2: ReadStats,
    post_local_r1: ReadStats,
    post_local_r2: ReadStats,
    pair_verdict: FilterResult,
    write: Option<(WriteRec, WriteRec)>,
}

/// Pure per-pair transform — same shape as [`process_se_record`] but for
/// a pair. Pair-level verdict resolves as `Pass` iff both mates pass.
/// Returns owned buffers so the worker can hand them back across the
/// rayon boundary.
#[allow(clippy::too_many_lines, clippy::needless_pass_by_value)]
fn process_pe_record(
    rec: OwnedPeRecord,
    cfg: FilterConfig,
    adapter: Option<&AdapterConfig>,
    polyg: Option<PolyGConfig>,
    umi: Option<UmiConfig>,
) -> ProcessedPe {
    let OwnedPeRecord { rec1, rec2 } = rec;
    let mut pre_local_r1 = ReadStats::default();
    let mut pre_local_r2 = ReadStats::default();
    pre_local_r1.observe(&rec1.seq, &rec1.qual);
    pre_local_r2.observe(&rec2.seq, &rec2.qual);

    let (id1_buf, id2_buf, off1, off2) = if let Some(u) = umi {
        let donor = if u.loc == UmiLoc::Read1 {
            &rec1.seq
        } else {
            &rec2.seq
        };
        match umi_extract(&rec1.id, donor, u) {
            Some((new_id, off)) => {
                let id2 = stamp_umi(&rec2.id, &new_id[rec1.id.len()..]);
                let (o1, o2) = if u.loc == UmiLoc::Read1 {
                    (off, 0)
                } else {
                    (0, off)
                };
                (new_id, id2, o1, o2)
            }
            None => {
                return ProcessedPe {
                    pre_local_r1,
                    pre_local_r2,
                    post_local_r1: ReadStats::default(),
                    post_local_r2: ReadStats::default(),
                    pair_verdict: FilterResult::TooShort,
                    write: None,
                };
            }
        }
    } else {
        (rec1.id.clone(), rec2.id.clone(), 0, 0)
    };

    let seq1_u = &rec1.seq[off1..];
    let q1_u = &rec1.qual[off1..];
    let seq2_u = &rec2.seq[off2..];
    let q2_u = &rec2.qual[off2..];

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

    let (post_local_r1, post_local_r2, write) = if matches!(pair_verdict, FilterResult::Pass) {
        let mut p1 = ReadStats::default();
        let mut p2 = ReadStats::default();
        p1.observe(seq1_t, q1_t);
        p2.observe(seq2_t, q2_t);
        (
            p1,
            p2,
            Some((
                (id1_buf, seq1_t.to_vec(), q1_t.to_vec()),
                (id2_buf, seq2_t.to_vec(), q2_t.to_vec()),
            )),
        )
    } else {
        (ReadStats::default(), ReadStats::default(), None)
    };

    ProcessedPe {
        pre_local_r1,
        pre_local_r2,
        post_local_r1,
        post_local_r2,
        pair_verdict,
        write,
    }
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
    let mut chunk: Vec<OwnedPeRecord> = Vec::with_capacity(CHUNK_RECORDS);

    let mut done = false;
    loop {
        chunk.clear();
        while chunk.len() < CHUNK_RECORDS {
            let r1 = r1_reader.next();
            let r2 = r2_reader.next();
            match (r1, r2) {
                (Some(rec1), Some(rec2)) => {
                    let rec1 = rec1.map_err(|e| parse_err("malformed R1 record", e))?;
                    let rec2 = rec2.map_err(|e| parse_err("malformed R2 record", e))?;
                    let q1 = rec1
                        .qual()
                        .ok_or_else(|| RsomicsError::InvalidInput("R1 missing quality".into()))?;
                    let q2 = rec2
                        .qual()
                        .ok_or_else(|| RsomicsError::InvalidInput("R2 missing quality".into()))?;
                    chunk.push(OwnedPeRecord {
                        rec1: OwnedSeRecord {
                            id: rec1.id().to_vec(),
                            seq: rec1.seq().into_owned(),
                            qual: q1.to_vec(),
                        },
                        rec2: OwnedSeRecord {
                            id: rec2.id().to_vec(),
                            seq: rec2.seq().into_owned(),
                            qual: q2.to_vec(),
                        },
                    });
                }
                // End-of-stream: stop filling the chunk but still process
                // whatever partial content it already has before exiting.
                (None, None) => {
                    done = true;
                    break;
                }
                (Some(Err(e)), _) => return Err(parse_err("malformed R1 record", e)),
                (_, Some(Err(e))) => return Err(parse_err("malformed R2 record", e)),
                (Some(_), None) | (None, Some(_)) => {
                    return Err(RsomicsError::InvalidInput(format!(
                        "paired-end inputs have different record counts: {} vs {}",
                        in1.display(),
                        in2.display(),
                    )));
                }
            }
        }
        if chunk.is_empty() {
            break;
        }

        let processed: Vec<ProcessedPe> = chunk
            .par_drain(..)
            .map(|rec| process_pe_record(rec, cfg, adapter, polyg, umi))
            .collect();

        for p in processed {
            pre_r1.merge(&p.pre_local_r1);
            pre_r2.merge(&p.pre_local_r2);
            // PE: each mate is one read in the filter counts (fastp compat).
            filtering.record(p.pair_verdict);
            filtering.record(p.pair_verdict);
            if let Some((a, b)) = p.write {
                post_r1.merge(&p.post_local_r1);
                post_r2.merge(&p.post_local_r2);
                write_record(&mut w1, &a.0, &a.1, &a.2).rs_context("writing R1 record")?;
                write_record(&mut w2, &b.0, &b.1, &b.2).rs_context("writing R2 record")?;
            }
        }
        if done {
            break;
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
