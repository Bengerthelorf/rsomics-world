// u64→f64 cast: lengths/counts fit in 52-bit mantissa for any real genome;
// cast only at final quartile/N50/percentage stage.
#![allow(clippy::cast_precision_loss)]

use std::path::Path;

use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};
use rsomics_seqstats::{LengthStats, classify, count_any_of};
use serde::Serialize;

pub use rsomics_seqstats::SeqType;

#[derive(Debug, Clone)]
pub struct Config {
    pub extended: bool,
    pub gap_letters: Vec<u8>,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            extended: false,
            gap_letters: b"- .".to_vec(),
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct FastaStats {
    pub file: String,
    pub format: &'static str,
    #[serde(rename = "type")]
    pub seq_type: SeqType,
    pub num_seqs: u64,
    pub sum_len: u64,
    pub min_len: u64,
    pub max_len: u64,
    pub avg_len: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub extended: Option<ExtendedStats>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ExtendedStats {
    #[serde(rename = "Q1")]
    pub q1: f64,
    #[serde(rename = "Q2")]
    pub q2: f64,
    #[serde(rename = "Q3")]
    pub q3: f64,
    pub sum_gap: u64,
    #[serde(rename = "N50")]
    pub n50: u64,
    // seqkit calls this N50_num; named L50 here. --tabular renders N50_num for compat.
    #[serde(rename = "L50")]
    pub l50: u64,
    #[serde(rename = "GC(%)")]
    pub gc_percent: f64,
    pub sum_n: u64,
}

const ALPHABET_GUESS_LIMIT: usize = 10_000;

#[allow(clippy::missing_errors_doc)]
pub fn compute_stats(path: &Path, cfg: &Config) -> Result<FastaStats> {
    let mut reader = parse_fastx_file(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("opening {}: {e}", path.display())))?;

    let mut lengths: Vec<u64> = Vec::new();
    let mut num_seqs: u64 = 0;
    let mut sum_len: u64 = 0;
    let mut min_len: u64 = u64::MAX;
    let mut max_len: u64 = 0;
    let mut sum_gap: u64 = 0;
    let mut sum_gc: u64 = 0;
    let mut sum_n_nuc: u64 = 0;
    let mut sum_n_prot: u64 = 0;

    let mut alphabet_sample: Vec<u8> = Vec::with_capacity(ALPHABET_GUESS_LIMIT);

    while let Some(record) = reader.next() {
        let rec = record
            .map_err(|e| RsomicsError::InvalidInput(format!("parsing {}: {e}", path.display())))?;
        let seq_cow = rec.seq();
        let seq: &[u8] = &seq_cow;
        let len = seq.len() as u64;

        num_seqs += 1;
        sum_len += len;
        if len < min_len {
            min_len = len;
        }
        if len > max_len {
            max_len = len;
        }
        if cfg.extended {
            lengths.push(len);
        }

        if alphabet_sample.len() < ALPHABET_GUESS_LIMIT {
            let take = (ALPHABET_GUESS_LIMIT - alphabet_sample.len()).min(seq.len());
            alphabet_sample.extend_from_slice(&seq[..take]);
        }

        if cfg.extended {
            sum_gap += count_any_of(seq, &cfg.gap_letters);
            sum_gc += count_any_of(seq, b"GCgc");
            sum_n_nuc += count_any_of(seq, b"Nn");
            sum_n_prot += count_any_of(seq, b"Xx");
        }
    }

    if num_seqs == 0 {
        return Err(RsomicsError::InvalidInput(format!(
            "{} contained no FASTA records",
            path.display()
        )));
    }

    let avg_len = sum_len as f64 / num_seqs as f64;
    let seq_type = classify(&alphabet_sample);

    let extended = if cfg.extended {
        let n_count = match seq_type {
            SeqType::Protein => sum_n_prot,
            _ => sum_n_nuc,
        };
        Some(extend(
            &mut lengths,
            sum_len,
            sum_gap,
            sum_gc,
            n_count,
            seq_type,
        ))
    } else {
        None
    };

    Ok(FastaStats {
        file: path.display().to_string(),
        format: "FASTA",
        seq_type,
        num_seqs,
        sum_len,
        min_len,
        max_len,
        avg_len,
        extended,
    })
}

fn extend(
    lengths: &mut Vec<u64>,
    sum_len: u64,
    sum_gap: u64,
    sum_gc: u64,
    sum_n: u64,
    seq_type: SeqType,
) -> ExtendedStats {
    let ls = LengthStats::new(std::mem::take(lengths));
    let (q1, q2, q3) = (ls.q1(), ls.q2(), ls.q3());
    let (n50, l50) = ls.n50_l50();
    let gc_percent = if matches!(seq_type, SeqType::Protein) || sum_len == 0 {
        0.0
    } else {
        sum_gc as f64 * 100.0 / sum_len as f64
    };
    ExtendedStats {
        q1,
        q2,
        q3,
        sum_gap,
        n50,
        l50,
        gc_percent,
        sum_n,
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    fn write_fa(s: &str) -> tempfile::NamedTempFile {
        use std::io::Write;
        let mut f = tempfile::Builder::new()
            .suffix(".fa")
            .tempfile()
            .expect("tempfile");
        f.write_all(s.as_bytes()).expect("write");
        f
    }

    #[test]
    fn empty_input_errors() {
        let f = write_fa("");
        let err = compute_stats(f.path(), &Config::default()).unwrap_err();
        assert!(matches!(err, RsomicsError::InvalidInput(_)));
    }

    #[test]
    fn basic_counts() {
        let f = write_fa(">a\nACGT\n>b\nACGTNN\n>c\nGGGG\n");
        let s = compute_stats(f.path(), &Config::default()).unwrap();
        assert_eq!(s.num_seqs, 3);
        assert_eq!(s.sum_len, 14);
        assert_eq!(s.min_len, 4);
        assert_eq!(s.max_len, 6);
        assert!((s.avg_len - 14.0 / 3.0).abs() < 1e-9);
        assert_eq!(s.seq_type, SeqType::Dna);
    }

    #[test]
    fn extended_gc_and_n() {
        let f = write_fa(">a\nACGT\n>b\nACGTNN\n>c\nGGGG\n");
        let cfg = Config {
            extended: true,
            ..Config::default()
        };
        let s = compute_stats(f.path(), &cfg).unwrap();
        let e = s.extended.expect("--all requested");
        assert_eq!(e.sum_gap, 0);
        assert_eq!(e.sum_n, 2);
        assert!((e.gc_percent - 8.0 / 14.0 * 100.0).abs() < 1e-9);
    }

    #[test]
    fn n50_three_contigs() {
        let f = write_fa(">a\nAAAA\n>b\nCCCCCC\n>c\nGGGGGGGG\n");
        let cfg = Config {
            extended: true,
            ..Config::default()
        };
        let s = compute_stats(f.path(), &cfg).unwrap();
        let e = s.extended.expect("--all");
        assert_eq!(e.n50, 6);
        assert_eq!(e.l50, 2);
    }

    #[test]
    fn rna_detection() {
        let f = write_fa(">r\nACGU\n");
        let s = compute_stats(f.path(), &Config::default()).unwrap();
        assert_eq!(s.seq_type, SeqType::Rna);
    }

    #[test]
    fn protein_detection() {
        let f = write_fa(">p\nMEEPSILQRT\n");
        let s = compute_stats(f.path(), &Config::default()).unwrap();
        assert_eq!(s.seq_type, SeqType::Protein);
    }

    #[test]
    fn gap_counting_uses_configured_letters() {
        let f = write_fa(">a\nA-C..G\n>b\nA C G\n");
        let cfg = Config {
            extended: true,
            ..Config::default()
        };
        let s = compute_stats(f.path(), &cfg).unwrap();
        let e = s.extended.expect("--all");
        assert_eq!(e.sum_gap, 5);
    }
}
