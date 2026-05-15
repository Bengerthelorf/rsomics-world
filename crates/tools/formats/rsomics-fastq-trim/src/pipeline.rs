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
    pub compression: i32,
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
            compression: 4,
        }
    }
}

/// Counters returned by [`Pipeline::run_se`] / [`Pipeline::run_pe`].
/// Serialised inside the `--json` envelope's `result`.
#[derive(Debug, Default, Clone, Serialize)]
pub struct TrimReport {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mode: Option<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub input_r1: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub input_r2: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output_r1: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output_r2: Option<String>,
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

impl std::ops::AddAssign<&TrimReport> for TrimReport {
    fn add_assign(&mut self, other: &TrimReport) {
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
        let mut writer = ChunkedWriter::create(output, self.cfg.compression)?;

        let mut report = TrimReport {
            mode: Some("SE"),
            input_r1: Some(input.display().to_string()),
            output_r1: Some(output.display().to_string()),
            ..TrimReport::default()
        };
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
                report += &p.delta;
                if let Some(t) = p.write {
                    writer.write_record(&t.id, t.seq_window(), t.qual_window())?;
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
        let mut w1 = ChunkedWriter::create(out1, self.cfg.compression)?;
        let mut w2 = ChunkedWriter::create(out2, self.cfg.compression)?;

        let mut report = TrimReport {
            mode: Some("PE"),
            input_r1: Some(in1.display().to_string()),
            input_r2: Some(in2.display().to_string()),
            output_r1: Some(out1.display().to_string()),
            output_r2: Some(out2.display().to_string()),
            ..TrimReport::default()
        };
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
                report += &p.delta;
                if let Some((t1, t2)) = p.write {
                    w1.write_record(&t1.id, t1.seq_window(), t1.qual_window())?;
                    w2.write_record(&t2.id, t2.seq_window(), t2.qual_window())?;
                }
            }
        }
        w1.finalize()?;
        w2.finalize()?;
        Ok(report)
    }
}

/// A trimmed record kept in its original-allocation `Vec`s with a live
/// `[start, end)` window. Skips the O(n) 5'-shift that `Vec::drain(..start)`
/// would force on every front-trimmed record.
struct TrimmedRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    start: usize,
    end: usize,
}

impl TrimmedRecord {
    fn seq_window(&self) -> &[u8] {
        &self.seq[self.start..self.end]
    }
    fn qual_window(&self) -> &[u8] {
        &self.qual[self.start..self.end]
    }
}

struct ProcessedSe {
    delta: TrimReport,
    write: Option<TrimmedRecord>,
}

struct ProcessedPe {
    delta: TrimReport,
    write: Option<(TrimmedRecord, TrimmedRecord)>,
}

#[allow(clippy::needless_pass_by_value)]
fn trim_se_record(rec: OwnedRecord, cfg: &PipelineConfig) -> ProcessedSe {
    let mut delta = TrimReport {
        reads_in: 1,
        bases_in: rec.seq.len() as u64,
        ..Default::default()
    };

    let original_len = rec.seq.len();

    // Operate on a sliding `[start, end)` window into rec.seq / rec.qual.
    // Avoids the O(n) `Vec::drain(..start)` shift on every fixed-trim run.
    let Some((start, mut end)) = apply_fixed(original_len, cfg.fixed1) else {
        delta.reads_too_short_after_trim = 1;
        return ProcessedSe { delta, write: None };
    };
    delta.fixed_trimmed_bases = (start + (original_len - end)) as u64;

    let window_len = end - start;
    let cut_adapter = cfg
        .adapter1
        .as_ref()
        .and_then(|a| find_adapter_3p(&rec.seq[start..end], a))
        .unwrap_or(window_len);
    if cut_adapter < window_len {
        delta.adapter_trimmed_reads = 1;
        delta.adapter_trimmed_bases = (window_len - cut_adapter) as u64;
        end = start + cut_adapter;
    }

    if let Some(pg) = cfg.poly_g
        && let Some(cut) = find_polyx_3p(&rec.seq[start..end], pg)
    {
        delta.poly_g_trimmed_reads = 1;
        delta.poly_g_trimmed_bases = ((end - start) - cut) as u64;
        end = start + cut;
    }

    if let Some(px) = cfg.poly_x
        && let Some(r) = find_dominant_polyx_3p(&rec.seq[start..end], px)
    {
        delta.poly_x_trimmed_reads = 1;
        delta.poly_x_trimmed_bases = ((end - start) - r.trim_at) as u64;
        end = start + r.trim_at;
    }

    if end - start < cfg.min_length_required {
        delta.reads_too_short_after_trim = 1;
        return ProcessedSe { delta, write: None };
    }

    delta.reads_out = 1;
    delta.bases_out = (end - start) as u64;

    ProcessedSe {
        delta,
        write: Some(TrimmedRecord {
            id: rec.id,
            seq: rec.seq,
            qual: rec.qual,
            start,
            end,
        }),
    }
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn trim_pe_pair(pair: OwnedPair, cfg: &PipelineConfig) -> ProcessedPe {
    let OwnedPair { r1, r2 } = pair;
    let orig1 = r1.seq.len();
    let orig2 = r2.seq.len();
    let mut delta = TrimReport {
        reads_in: 2,
        bases_in: (orig1 + orig2) as u64,
        ..Default::default()
    };

    let Some((s1, mut e1)) = apply_fixed(orig1, cfg.fixed1) else {
        delta.reads_too_short_after_trim = 2;
        return ProcessedPe { delta, write: None };
    };
    let Some((s2, mut e2)) = apply_fixed(orig2, cfg.fixed2) else {
        delta.reads_too_short_after_trim = 2;
        return ProcessedPe { delta, write: None };
    };
    delta.fixed_trimmed_bases = (s1 + (orig1 - e1) + s2 + (orig2 - e2)) as u64;

    // PolyG must run BEFORE overlap analysis (fastp PE order).
    if let Some(pg) = cfg.poly_g {
        if let Some(cut) = find_polyx_3p(&r1.seq[s1..e1], pg) {
            delta.poly_g_trimmed_reads += 1;
            delta.poly_g_trimmed_bases += ((e1 - s1) - cut) as u64;
            e1 = s1 + cut;
        }
        if let Some(cut) = find_polyx_3p(&r2.seq[s2..e2], pg) {
            delta.poly_g_trimmed_reads += 1;
            delta.poly_g_trimmed_bases += ((e2 - s2) - cut) as u64;
            e2 = s2 + cut;
        }
    }

    let mut overlap_fired = false;
    if let Some(ov_cfg) = cfg.overlap {
        let r2_rc = reverse_complement(&r2.seq[s2..e2]);
        let ov: OverlapResult = overlap_analyze(&r1.seq[s1..e1], &r2_rc, ov_cfg);
        if let Some((new1, new2)) = overlap_trim_lengths(ov, e1 - s1, e2 - s2, s1, s2) {
            delta.overlap_trimmed_pairs = 1;
            delta.overlap_trimmed_bases += (((e1 - s1) - new1) + ((e2 - s2) - new2)) as u64;
            e1 = s1 + new1;
            e2 = s2 + new2;
            overlap_fired = true;
        }
    }

    if !overlap_fired {
        if let Some(a1) = cfg.adapter1.as_ref()
            && let Some(cut) = find_adapter_3p(&r1.seq[s1..e1], a1)
        {
            delta.adapter_trimmed_reads += 1;
            delta.adapter_trimmed_bases += ((e1 - s1) - cut) as u64;
            e1 = s1 + cut;
        }
        if let Some(a2) = cfg.adapter2.as_ref()
            && let Some(cut) = find_adapter_3p(&r2.seq[s2..e2], a2)
        {
            delta.adapter_trimmed_reads += 1;
            delta.adapter_trimmed_bases += ((e2 - s2) - cut) as u64;
            e2 = s2 + cut;
        }
    }

    if let Some(px) = cfg.poly_x {
        if let Some(r) = find_dominant_polyx_3p(&r1.seq[s1..e1], px) {
            delta.poly_x_trimmed_reads += 1;
            delta.poly_x_trimmed_bases += ((e1 - s1) - r.trim_at) as u64;
            e1 = s1 + r.trim_at;
        }
        if let Some(r) = find_dominant_polyx_3p(&r2.seq[s2..e2], px) {
            delta.poly_x_trimmed_reads += 1;
            delta.poly_x_trimmed_bases += ((e2 - s2) - r.trim_at) as u64;
            e2 = s2 + r.trim_at;
        }
    }

    if (e1 - s1) < cfg.min_length_required || (e2 - s2) < cfg.min_length_required {
        delta.reads_too_short_after_trim = 2;
        return ProcessedPe { delta, write: None };
    }

    delta.reads_out = 2;
    delta.bases_out = ((e1 - s1) + (e2 - s2)) as u64;

    ProcessedPe {
        delta,
        write: Some((
            TrimmedRecord {
                id: r1.id,
                seq: r1.seq,
                qual: r1.qual,
                start: s1,
                end: e1,
            },
            TrimmedRecord {
                id: r2.id,
                seq: r2.seq,
                qual: r2.qual,
                start: s2,
                end: e2,
            },
        )),
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
