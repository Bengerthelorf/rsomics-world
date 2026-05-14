//! Rayon-chunked scatter/gather over a FASTQ input.
//!
//! Per-record work order (matches fastp `seprocessor.cpp` for SE,
//! `peprocessor.cpp` for PE, with the per-function partition trimming
//! the quality / UMI / stats stages that live in sibling crates):
//!
//! - SE: fixed → adapter → polyG → polyX → emit
//! - PE: fixed → polyG → overlap → (fallback static adapter) → polyX → emit
//!
//! The only quality-adjacent check kept in this crate is the **post-trim
//! min-length gate** — a zero-length read after trim must be discarded
//! since no downstream tool will emit it. Full quality / N-content /
//! sliding-window filtering lives in `rsomics-fastq-quality`.

use std::path::Path;

use needletail::parse_fastx_file;
use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::adapter::{AdapterConfig, find_adapter_3p};
use crate::fixed::{FixedTrimConfig, apply_fixed};
use crate::overlap::{
    OverlapConfig, OverlapResult, analyze as overlap_analyze, reverse_complement,
    trim_lengths as overlap_trim_lengths,
};
use crate::parallel_gz::ChunkedWriter;
use crate::polyx::{PolyXConfig, find_dominant_polyx_3p, find_polyx_3p};

/// Chunk size for the parallel scatter/gather. Larger amortises rayon
/// dispatch overhead; smaller reduces memory peak. 8192 records ≈ 12 MB
/// of sequence per chunk for typical 150 bp reads — comfortable on any
/// modern machine and significantly fewer dispatches per file.
const CHUNK_RECORDS: usize = 8192;

/// One FASTQ record decoupled from needletail's borrowed buffers.
struct OwnedRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

struct OwnedPair {
    r1: OwnedRecord,
    r2: OwnedRecord,
}

/// Trim configuration assembled from CLI args, shared across rayon
/// workers via reference.
#[derive(Debug, Clone)]
pub struct PipelineConfig {
    pub fixed1: FixedTrimConfig,
    pub fixed2: FixedTrimConfig,
    pub adapter1: Option<AdapterConfig>,
    pub adapter2: Option<AdapterConfig>,
    pub poly_g: Option<PolyXConfig>,
    pub poly_x: Option<PolyXConfig>,
    pub overlap: Option<OverlapConfig>,
    pub min_length_required: usize,
}

impl Default for PipelineConfig {
    fn default() -> Self {
        Self {
            fixed1: FixedTrimConfig::default(),
            fixed2: FixedTrimConfig::default(),
            adapter1: None,
            adapter2: None,
            poly_g: None,
            poly_x: None,
            overlap: None,
            min_length_required: 15,
        }
    }
}

/// Counters returned by [`Pipeline::run_se`] / [`Pipeline::run_pe`].
/// Serialised inside the `--json` envelope's `result`.
#[derive(Debug, Default, Clone, Serialize)]
pub struct TrimReport {
    pub reads_in: u64,
    pub reads_out: u64,
    pub bases_in: u64,
    pub bases_out: u64,
    pub adapter_trimmed_reads: u64,
    pub adapter_trimmed_bases: u64,
    pub poly_g_trimmed_reads: u64,
    pub poly_g_trimmed_bases: u64,
    pub poly_x_trimmed_reads: u64,
    pub poly_x_trimmed_bases: u64,
    pub fixed_trimmed_bases: u64,
    pub overlap_trimmed_pairs: u64,
    pub overlap_trimmed_bases: u64,
    pub reads_too_short_after_trim: u64,
}

impl TrimReport {
    fn merge(&mut self, other: &Self) {
        self.reads_in += other.reads_in;
        self.reads_out += other.reads_out;
        self.bases_in += other.bases_in;
        self.bases_out += other.bases_out;
        self.adapter_trimmed_reads += other.adapter_trimmed_reads;
        self.adapter_trimmed_bases += other.adapter_trimmed_bases;
        self.poly_g_trimmed_reads += other.poly_g_trimmed_reads;
        self.poly_g_trimmed_bases += other.poly_g_trimmed_bases;
        self.poly_x_trimmed_reads += other.poly_x_trimmed_reads;
        self.poly_x_trimmed_bases += other.poly_x_trimmed_bases;
        self.fixed_trimmed_bases += other.fixed_trimmed_bases;
        self.overlap_trimmed_pairs += other.overlap_trimmed_pairs;
        self.overlap_trimmed_bases += other.overlap_trimmed_bases;
        self.reads_too_short_after_trim += other.reads_too_short_after_trim;
    }
}

/// Pipeline shape: owns the `PipelineConfig`, exposes `run_se` and
/// `run_pe` as the two operational entry points.
pub struct Pipeline<'cfg> {
    pub cfg: &'cfg PipelineConfig,
}

impl<'cfg> Pipeline<'cfg> {
    #[must_use]
    pub fn new(cfg: &'cfg PipelineConfig) -> Self {
        Self { cfg }
    }

    /// Stream a single-end FASTQ through the trim pipeline.
    ///
    /// # Errors
    ///
    /// `InvalidInput` if input parsing fails; `Io` if output write fails.
    pub fn run_se(&self, input: &Path, output: &Path) -> Result<TrimReport> {
        let mut reader = parse_fastx_file(input)
            .map_err(|e| parse_err(&format!("opening input {}", input.display()), e))?;
        let mut writer = ChunkedWriter::create(output)?;

        let mut report = TrimReport::default();
        let mut chunk: Vec<OwnedRecord> = Vec::with_capacity(CHUNK_RECORDS);

        loop {
            chunk.clear();
            while chunk.len() < CHUNK_RECORDS {
                let Some(r) = reader.next() else { break };
                let rec = r.map_err(|e| parse_err("malformed FASTQ record", e))?;
                let qual = rec.qual().ok_or_else(|| {
                    RsomicsError::InvalidInput("FASTQ record missing quality line".into())
                })?;
                chunk.push(OwnedRecord {
                    id: rec.id().to_vec(),
                    seq: rec.seq().into_owned(),
                    qual: qual.to_vec(),
                });
            }
            if chunk.is_empty() {
                break;
            }

            let processed: Vec<ProcessedSe> = chunk
                .par_drain(..)
                .map(|rec| trim_se_record(rec, self.cfg))
                .collect();

            for p in processed {
                report.merge(&p.delta);
                if let Some((id, seq, qual)) = p.write {
                    writer.write_record(&id, &seq, &qual)?;
                }
            }
        }
        writer.finalize()?;
        Ok(report)
    }

    /// Stream a paired-end FASTQ through the trim pipeline.
    ///
    /// # Errors
    ///
    /// `InvalidInput` if either input parses incorrectly or the pair
    /// counts diverge; `Io` if writes fail.
    pub fn run_pe(&self, in1: &Path, in2: &Path, out1: &Path, out2: &Path) -> Result<TrimReport> {
        let mut r1 = parse_fastx_file(in1)
            .map_err(|e| parse_err(&format!("opening input {}", in1.display()), e))?;
        let mut r2 = parse_fastx_file(in2)
            .map_err(|e| parse_err(&format!("opening input {}", in2.display()), e))?;
        let mut w1 = ChunkedWriter::create(out1)?;
        let mut w2 = ChunkedWriter::create(out2)?;

        let mut report = TrimReport::default();
        let mut chunk: Vec<OwnedPair> = Vec::with_capacity(CHUNK_RECORDS);

        let mut done = false;
        while !done {
            chunk.clear();
            while chunk.len() < CHUNK_RECORDS {
                let (a, b) = (r1.next(), r2.next());
                match (a, b) {
                    (Some(ra), Some(rb)) => {
                        let rec1 = own_record(ra)?;
                        let rec2 = own_record(rb)?;
                        chunk.push(OwnedPair { r1: rec1, r2: rec2 });
                    }
                    (None, None) => {
                        done = true;
                        break;
                    }
                    _ => {
                        return Err(RsomicsError::InvalidInput(
                            "PE input record counts diverge".into(),
                        ));
                    }
                }
            }
            if chunk.is_empty() {
                break;
            }

            let processed: Vec<ProcessedPe> = chunk
                .par_drain(..)
                .map(|pair| trim_pe_pair(pair, self.cfg))
                .collect();

            for p in processed {
                report.merge(&p.delta);
                if let Some((rec1, rec2)) = p.write {
                    w1.write_record(&rec1.0, &rec1.1, &rec1.2)?;
                    w2.write_record(&rec2.0, &rec2.1, &rec2.2)?;
                }
            }
        }
        w1.finalize()?;
        w2.finalize()?;
        Ok(report)
    }
}

type WriteRec = (Vec<u8>, Vec<u8>, Vec<u8>);

struct ProcessedSe {
    delta: TrimReport,
    write: Option<WriteRec>,
}

struct ProcessedPe {
    delta: TrimReport,
    write: Option<(WriteRec, WriteRec)>,
}

#[allow(clippy::needless_pass_by_value)]
fn trim_se_record(rec: OwnedRecord, cfg: &PipelineConfig) -> ProcessedSe {
    let mut delta = TrimReport {
        reads_in: 1,
        bases_in: rec.seq.len() as u64,
        ..Default::default()
    };

    let original_len = rec.seq.len();

    // 1. Fixed-length trim.
    let Some((start, end)) = apply_fixed(rec.seq.len(), cfg.fixed1) else {
        delta.reads_too_short_after_trim = 1;
        return ProcessedSe { delta, write: None };
    };
    let fixed_bases = (start + (original_len - end)) as u64;
    delta.fixed_trimmed_bases = fixed_bases;

    let seq = &rec.seq[start..end];
    let qual = &rec.qual[start..end];

    // 2. 3' adapter trim (static sequence, if configured).
    let after_adapter = cfg
        .adapter1
        .as_ref()
        .and_then(|a| find_adapter_3p(seq, a))
        .unwrap_or(seq.len());
    if after_adapter < seq.len() {
        delta.adapter_trimmed_reads = 1;
        delta.adapter_trimmed_bases = (seq.len() - after_adapter) as u64;
    }

    // 3. PolyG trim.
    let after_polyg = cfg
        .poly_g
        .and_then(|pg| find_polyx_3p(&seq[..after_adapter], pg))
        .unwrap_or(after_adapter);
    if after_polyg < after_adapter {
        delta.poly_g_trimmed_reads = 1;
        delta.poly_g_trimmed_bases = (after_adapter - after_polyg) as u64;
    }

    // 4. PolyX trim.
    let after_px = cfg
        .poly_x
        .and_then(|px| find_dominant_polyx_3p(&seq[..after_polyg], px))
        .map_or(after_polyg, |r| r.trim_at);
    if after_px < after_polyg {
        delta.poly_x_trimmed_reads = 1;
        delta.poly_x_trimmed_bases = (after_polyg - after_px) as u64;
    }

    let trim_at = after_px;
    if trim_at < cfg.min_length_required {
        delta.reads_too_short_after_trim = 1;
        return ProcessedSe { delta, write: None };
    }

    let seq_out = seq[..trim_at].to_vec();
    let qual_out = qual[..trim_at].to_vec();
    delta.reads_out = 1;
    delta.bases_out = seq_out.len() as u64;

    ProcessedSe {
        delta,
        write: Some((rec.id, seq_out, qual_out)),
    }
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn trim_pe_pair(pair: OwnedPair, cfg: &PipelineConfig) -> ProcessedPe {
    let OwnedPair { r1, r2 } = pair;
    let mut delta = TrimReport {
        reads_in: 2,
        bases_in: (r1.seq.len() + r2.seq.len()) as u64,
        ..Default::default()
    };

    // 1. Fixed-length trim for each mate independently.
    let Some((s1, e1)) = apply_fixed(r1.seq.len(), cfg.fixed1) else {
        delta.reads_too_short_after_trim = 2;
        return ProcessedPe { delta, write: None };
    };
    let Some((s2, e2)) = apply_fixed(r2.seq.len(), cfg.fixed2) else {
        delta.reads_too_short_after_trim = 2;
        return ProcessedPe { delta, write: None };
    };
    delta.fixed_trimmed_bases = (s1 + (r1.seq.len() - e1) + s2 + (r2.seq.len() - e2)) as u64;

    let mut seq1: Vec<u8> = r1.seq[s1..e1].to_vec();
    let mut qual1: Vec<u8> = r1.qual[s1..e1].to_vec();
    let mut seq2: Vec<u8> = r2.seq[s2..e2].to_vec();
    let mut qual2: Vec<u8> = r2.qual[s2..e2].to_vec();

    // 2. PolyG trim per mate (fastp does this BEFORE overlap analysis).
    if let Some(pg) = cfg.poly_g {
        if let Some(cut) = find_polyx_3p(&seq1, pg) {
            delta.poly_g_trimmed_reads += 1;
            delta.poly_g_trimmed_bases += (seq1.len() - cut) as u64;
            seq1.truncate(cut);
            qual1.truncate(cut);
        }
        if let Some(cut) = find_polyx_3p(&seq2, pg) {
            delta.poly_g_trimmed_reads += 1;
            delta.poly_g_trimmed_bases += (seq2.len() - cut) as u64;
            seq2.truncate(cut);
            qual2.truncate(cut);
        }
    }

    // 3. PE overlap analysis. RC of R2 is built once per pair.
    let mut overlap_fired = false;
    if let Some(ov_cfg) = cfg.overlap {
        let r2_rc = reverse_complement(&seq2);
        let ov: OverlapResult = overlap_analyze(&seq1, &r2_rc, ov_cfg);
        if let Some((new1, new2)) = overlap_trim_lengths(ov, seq1.len(), seq2.len(), s1, s2) {
            delta.overlap_trimmed_pairs = 1;
            delta.overlap_trimmed_bases += ((seq1.len() - new1) + (seq2.len() - new2)) as u64;
            seq1.truncate(new1);
            qual1.truncate(new1);
            seq2.truncate(new2);
            qual2.truncate(new2);
            overlap_fired = true;
        }
    }

    // 4. Static-sequence adapter trim (fallback when overlap didn't fire).
    if !overlap_fired {
        if let Some(a1) = cfg.adapter1.as_ref()
            && let Some(cut) = find_adapter_3p(&seq1, a1)
        {
            delta.adapter_trimmed_reads += 1;
            delta.adapter_trimmed_bases += (seq1.len() - cut) as u64;
            seq1.truncate(cut);
            qual1.truncate(cut);
        }
        if let Some(a2) = cfg.adapter2.as_ref()
            && let Some(cut) = find_adapter_3p(&seq2, a2)
        {
            delta.adapter_trimmed_reads += 1;
            delta.adapter_trimmed_bases += (seq2.len() - cut) as u64;
            seq2.truncate(cut);
            qual2.truncate(cut);
        }
    }

    if let Some(px) = cfg.poly_x {
        if let Some(r) = find_dominant_polyx_3p(&seq1, px) {
            delta.poly_x_trimmed_reads += 1;
            delta.poly_x_trimmed_bases += (seq1.len() - r.trim_at) as u64;
            seq1.truncate(r.trim_at);
            qual1.truncate(r.trim_at);
        }
        if let Some(r) = find_dominant_polyx_3p(&seq2, px) {
            delta.poly_x_trimmed_reads += 1;
            delta.poly_x_trimmed_bases += (seq2.len() - r.trim_at) as u64;
            seq2.truncate(r.trim_at);
            qual2.truncate(r.trim_at);
        }
    }

    if seq1.len() < cfg.min_length_required || seq2.len() < cfg.min_length_required {
        delta.reads_too_short_after_trim = 2;
        return ProcessedPe { delta, write: None };
    }

    delta.reads_out = 2;
    delta.bases_out = (seq1.len() + seq2.len()) as u64;

    ProcessedPe {
        delta,
        write: Some(((r1.id, seq1, qual1), (r2.id, seq2, qual2))),
    }
}

fn own_record(
    r: std::result::Result<needletail::parser::SequenceRecord, needletail::errors::ParseError>,
) -> Result<OwnedRecord> {
    let rec = r.map_err(|e| parse_err("malformed FASTQ record", e))?;
    let qual = rec
        .qual()
        .ok_or_else(|| RsomicsError::InvalidInput("FASTQ record missing quality line".into()))?;
    Ok(OwnedRecord {
        id: rec.id().to_vec(),
        seq: rec.seq().into_owned(),
        qual: qual.to_vec(),
    })
}

fn parse_err(prefix: &str, e: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!("{prefix}: {e}"))
}
