use std::path::Path;

use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::{OwnedRecord, open_fastq};
use serde::Serialize;

use rsomics_fqgz::ChunkedWriter;

const CHUNK_RECORDS: usize = 8192; // ≈12 MB per chunk at 150 bp; amortises rayon dispatch

// fastp check order (src/filter.cpp Filter::passFilter): low-qual-percent → n-base → length.
// Order determines which failure counter increments when a read fails multiple criteria.
pub enum FilterOutcome {
    Pass,
    FailLength,
    FailLowQual,
    FailNBase,
}

// fastp src/filter.cpp + src/options.h (MIT; reading permitted).
// qual_threshold_ascii is u16 to prevent overflow: Phred 223 + 33 = 256 wraps in u8,
// emitting silent wrong output in release. The widened comparison is zero-cost.
#[derive(Debug, Clone, Copy)]
pub struct FilterConfig {
    // ASCII threshold: a quality byte strictly less than this is low-quality.
    pub qual_threshold_ascii: u16,
    pub unqualified_percent_limit: u8,
    pub n_base_limit: usize,
    pub required_length: usize,
    // 0 = no upper bound.
    pub max_length: usize,
    pub quality_enabled: bool,
    pub length_enabled: bool,
}

impl FilterConfig {
    #[must_use]
    pub fn check(&self, seq: &[u8], qual: &[u8]) -> FilterOutcome {
        let len = seq.len();

        if self.quality_enabled {
            let mut low_qual: usize = 0;
            let mut n_count: usize = 0;
            for (&b, &q) in seq.iter().zip(qual.iter()) {
                if u16::from(q) < self.qual_threshold_ascii {
                    low_qual += 1;
                }
                if b == b'N' || b == b'n' {
                    n_count += 1;
                }
            }
            // Integer-only rewrite of fastp's `lowQualNum > (unqualifiedPercentLimit * rlen / 100.0)`.
            // Equivalent because reads fit comfortably in usize; no fp precision concern.
            if low_qual * 100 > usize::from(self.unqualified_percent_limit) * len {
                return FilterOutcome::FailLowQual;
            }
            if n_count > self.n_base_limit {
                return FilterOutcome::FailNBase;
            }
        }

        if self.length_enabled {
            if len < self.required_length {
                return FilterOutcome::FailLength;
            }
            if self.max_length > 0 && len > self.max_length {
                return FilterOutcome::FailLength;
            }
        }

        FilterOutcome::Pass
    }
}

#[derive(Debug, Default, Clone, Serialize)]
pub struct FilterReport {
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
    pub reads_failed_quality: u64,
    pub reads_failed_length: u64,
    pub reads_failed_n_bases: u64,
}

impl std::ops::AddAssign<&FilterReport> for FilterReport {
    fn add_assign(&mut self, other: &FilterReport) {
        self.reads_in += other.reads_in;
        self.reads_out += other.reads_out;
        self.bases_in += other.bases_in;
        self.bases_out += other.bases_out;
        self.reads_failed_quality += other.reads_failed_quality;
        self.reads_failed_length += other.reads_failed_length;
        self.reads_failed_n_bases += other.reads_failed_n_bases;
    }
}

pub struct Pipeline<'cfg> {
    pub cfg: &'cfg FilterConfig,
    pub compression: i32,
}

impl<'cfg> Pipeline<'cfg> {
    #[must_use]
    pub fn new(cfg: &'cfg FilterConfig, compression: i32) -> Self {
        Self { cfg, compression }
    }

    #[allow(clippy::missing_errors_doc)]
    pub fn run_se(&self, input: &Path, output: &Path) -> Result<FilterReport> {
        let mut reader = open_fastq(input)?;
        let mut writer = ChunkedWriter::create(output, self.compression)?;

        let mut report = FilterReport {
            mode: Some("SE"),
            input_r1: Some(input.display().to_string()),
            output_r1: Some(output.display().to_string()),
            ..FilterReport::default()
        };
        let mut chunk: Vec<OwnedRecord> = Vec::with_capacity(CHUNK_RECORDS);

        loop {
            chunk.clear();
            while chunk.len() < CHUNK_RECORDS {
                let Some(r) = reader.next() else { break };
                chunk.push(r?);
            }
            if chunk.is_empty() {
                break;
            }

            let processed: Vec<ProcessedSe> = chunk
                .par_drain(..)
                .map(|rec| filter_se_record(rec, self.cfg))
                .collect();

            for p in processed {
                report += &p.delta;
                if let Some(rec) = p.write {
                    writer.write_record(&rec.id, &rec.seq, &rec.qual)?;
                }
            }
        }
        writer.finalize()?;
        Ok(report)
    }

    #[allow(clippy::missing_errors_doc)]
    pub fn run_pe(&self, in1: &Path, in2: &Path, out1: &Path, out2: &Path) -> Result<FilterReport> {
        let mut r1 = open_fastq(in1)?;
        let mut r2 = open_fastq(in2)?;
        let mut w1 = ChunkedWriter::create(out1, self.compression)?;
        let mut w2 = ChunkedWriter::create(out2, self.compression)?;

        let mut report = FilterReport {
            mode: Some("PE"),
            input_r1: Some(in1.display().to_string()),
            input_r2: Some(in2.display().to_string()),
            output_r1: Some(out1.display().to_string()),
            output_r2: Some(out2.display().to_string()),
            ..FilterReport::default()
        };
        let mut chunk: Vec<OwnedPair> = Vec::with_capacity(CHUNK_RECORDS);

        let mut done = false;
        while !done {
            chunk.clear();
            while chunk.len() < CHUNK_RECORDS {
                let (a, b) = (r1.next(), r2.next());
                match (a, b) {
                    (Some(ra), Some(rb)) => {
                        chunk.push(OwnedPair { r1: ra?, r2: rb? });
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
                .map(|pair| filter_pe_pair(pair, self.cfg))
                .collect();

            for p in processed {
                report += &p.delta;
                if let Some((rec1, rec2)) = p.write {
                    w1.write_record(&rec1.id, &rec1.seq, &rec1.qual)?;
                    w2.write_record(&rec2.id, &rec2.seq, &rec2.qual)?;
                }
            }
        }
        w1.finalize()?;
        w2.finalize()?;
        Ok(report)
    }
}

struct OwnedPair {
    r1: OwnedRecord,
    r2: OwnedRecord,
}

struct ProcessedSe {
    delta: FilterReport,
    write: Option<OwnedRecord>,
}

struct ProcessedPe {
    delta: FilterReport,
    write: Option<(OwnedRecord, OwnedRecord)>,
}

#[allow(clippy::needless_pass_by_value)]
fn filter_se_record(rec: OwnedRecord, cfg: &FilterConfig) -> ProcessedSe {
    let mut delta = FilterReport {
        reads_in: 1,
        bases_in: rec.seq.len() as u64,
        ..Default::default()
    };

    match cfg.check(&rec.seq, &rec.qual) {
        FilterOutcome::Pass => {
            delta.reads_out = 1;
            delta.bases_out = rec.seq.len() as u64;
            ProcessedSe {
                delta,
                write: Some(rec),
            }
        }
        FilterOutcome::FailLength => {
            delta.reads_failed_length = 1;
            ProcessedSe { delta, write: None }
        }
        FilterOutcome::FailLowQual => {
            delta.reads_failed_quality = 1;
            ProcessedSe { delta, write: None }
        }
        FilterOutcome::FailNBase => {
            delta.reads_failed_n_bases = 1;
            ProcessedSe { delta, write: None }
        }
    }
}

#[allow(clippy::needless_pass_by_value)]
fn filter_pe_pair(pair: OwnedPair, cfg: &FilterConfig) -> ProcessedPe {
    let OwnedPair { r1, r2 } = pair;
    let mut delta = FilterReport {
        reads_in: 2,
        bases_in: (r1.seq.len() + r2.seq.len()) as u64,
        ..Default::default()
    };

    let outcome1 = cfg.check(&r1.seq, &r1.qual);
    let outcome2 = cfg.check(&r2.seq, &r2.qual);

    let both_pass = matches!(
        (&outcome1, &outcome2),
        (FilterOutcome::Pass, FilterOutcome::Pass)
    );

    for outcome in [outcome1, outcome2] {
        match outcome {
            FilterOutcome::Pass => {}
            FilterOutcome::FailLength => delta.reads_failed_length += 1,
            FilterOutcome::FailLowQual => delta.reads_failed_quality += 1,
            FilterOutcome::FailNBase => delta.reads_failed_n_bases += 1,
        }
    }

    if both_pass {
        delta.reads_out = 2;
        delta.bases_out = (r1.seq.len() + r2.seq.len()) as u64;
        ProcessedPe {
            delta,
            write: Some((r1, r2)),
        }
    } else {
        ProcessedPe { delta, write: None }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config_both_enabled() -> FilterConfig {
        FilterConfig {
            // Q15 threshold: ASCII 48 ('0' in Phred+33)
            qual_threshold_ascii: 48,
            unqualified_percent_limit: 40,
            n_base_limit: 5,
            required_length: 20,
            max_length: 0,
            quality_enabled: true,
            length_enabled: true,
        }
    }

    // A read that is both too short AND has too many low-qual bases must yield
    // FailLowQual, not FailLength — matching fastp Filter::passFilter order:
    // low-qual-percent check precedes length check (fastp src/filter.cpp).
    #[test]
    fn quality_check_precedes_length_check() {
        let cfg = config_both_enabled();

        let seq = b"AAAAAAAAAA";
        let qual = b"!!!!!!!!!!"; // ASCII 33 = Q0, below threshold 48
        assert!(
            matches!(cfg.check(seq, qual), FilterOutcome::FailLowQual),
            "short + low-qual read must fail as FailLowQual (fastp quality check precedes length check)"
        );
    }

    #[test]
    fn short_read_good_quality_fails_length() {
        let cfg = config_both_enabled();

        let seq = b"AAAAAAAAAA";
        let qual = b"IIIIIIIIII";
        assert!(
            matches!(cfg.check(seq, qual), FilterOutcome::FailLength),
            "short read with good quality must fail as FailLength"
        );
    }
}
