// u64→f64 cast: lengths/counts fit in 52-bit mantissa for any real run;
// cast only at the final quartile/N50/percentage/AvgQual stage.
#![allow(clippy::cast_precision_loss)]

use std::path::Path;
use std::sync::LazyLock;

use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};
use rsomics_seqstats::{DEFAULT_ALPHABET_GUESS_LEN, LengthStats, classify, count_any_of};
use serde::Serialize;

pub use rsomics_seqstats::SeqType;

#[derive(Debug, Clone)]
pub struct Config {
    pub extended: bool,
    pub gap_letters: Vec<u8>,
    pub qual_offset: u8,
    pub basename: bool,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            extended: false,
            gap_letters: b"- .".to_vec(),
            qual_offset: 33,
            basename: false,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct FastqStats {
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
    #[serde(rename = "Q20(%)")]
    pub q20_percent: f64,
    #[serde(rename = "Q30(%)")]
    pub q30_percent: f64,
    #[serde(rename = "AvgQual")]
    pub avg_qual: f64,
    #[serde(rename = "GC(%)")]
    pub gc_percent: f64,
    pub sum_n: u64,
}

// seqkit rounds to 2 decimals first, then prints %.0f — two-step matches boundary values.
fn round2(x: f64) -> f64 {
    (x * 100.0).round() / 100.0
}

// 10^(-q/10) lookup by phred q; avoids per-base transcendental call — same approach as seqkit's qual_map.
static ERR_PROB: LazyLock<[f64; 256]> = LazyLock::new(|| {
    let mut t = [0.0f64; 256];
    for (q, e) in t.iter_mut().enumerate() {
        *e = 10f64.powf(-(q as f64) / 10.0);
    }
    t
});

#[allow(clippy::missing_errors_doc)]
pub fn compute_stats(path: &Path, cfg: &Config) -> Result<FastqStats> {
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
    let mut qc = QualCounts {
        q20: 0,
        q30: 0,
        err_sum: 0.0,
    };

    // seqkit guesses `type` from the first record only (its FASTX reader
    // sets the alphabet once); classify it inline so no per-record sample
    // is accumulated.
    let mut seq_type: Option<SeqType> = None;

    while let Some(record) = reader.next() {
        let rec = record
            .map_err(|e| RsomicsError::InvalidInput(format!("parsing {}: {e}", path.display())))?;
        let seq_cow = rec.seq();
        let seq: &[u8] = &seq_cow;
        let len = seq.len() as u64;

        let qual = rec.qual().ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "{}: record {} has no quality line — not a FASTQ",
                path.display(),
                num_seqs + 1
            ))
        })?;

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

        if seq_type.is_none() {
            let take = seq.len().min(DEFAULT_ALPHABET_GUESS_LEN);
            seq_type = Some(classify(&seq[..take]));
        }

        accumulate_qual(qual, cfg.qual_offset, &mut qc, path)?;

        if cfg.extended {
            sum_gap += count_any_of(seq, &cfg.gap_letters);
            sum_gc += count_any_of(seq, b"GCgc");
            sum_n_nuc += count_any_of(seq, b"Nn");
            sum_n_prot += count_any_of(seq, b"Xx");
        }
    }

    // `seq_type` is `None` exactly when there were no records — the
    // first record always sets it (even an empty read ⇒ Unlimit).
    let Some(seq_type) = seq_type else {
        return Err(RsomicsError::InvalidInput(format!(
            "{} contained no FASTQ records",
            path.display()
        )));
    };

    let avg_len = sum_len as f64 / num_seqs as f64;

    let extended = cfg.extended.then(|| {
        let sum_n = if matches!(seq_type, SeqType::Protein) {
            sum_n_prot
        } else {
            sum_n_nuc
        };
        extend(&mut lengths, seq_type, sum_len, sum_gap, sum_gc, sum_n, &qc)
    });

    let file = if cfg.basename {
        path.file_name().map_or_else(
            || path.display().to_string(),
            |n| n.to_string_lossy().into(),
        )
    } else {
        path.display().to_string()
    };

    Ok(FastqStats {
        file,
        format: "FASTQ",
        seq_type,
        num_seqs,
        sum_len,
        min_len,
        max_len,
        avg_len,
        extended,
    })
}

struct QualCounts {
    q20: u64,
    q30: u64,
    err_sum: f64,
}

fn accumulate_qual(qual: &[u8], offset: u8, qc: &mut QualCounts, path: &Path) -> Result<()> {
    for &qb in qual {
        let Some(q) = qb.checked_sub(offset) else {
            return Err(RsomicsError::InvalidInput(format!(
                "{}: quality byte {qb} below offset {offset} (wrong --fq-encoding?)",
                path.display()
            )));
        };
        if q >= 20 {
            qc.q20 += 1;
            if q >= 30 {
                qc.q30 += 1;
            }
        }
        qc.err_sum += ERR_PROB[q as usize];
    }
    Ok(())
}

fn extend(
    lengths: &mut Vec<u64>,
    seq_type: SeqType,
    sum_len: u64,
    sum_gap: u64,
    sum_gc: u64,
    sum_n: u64,
    qual: &QualCounts,
) -> ExtendedStats {
    let ls = LengthStats::new(std::mem::take(lengths));
    let (n50, l50) = ls.n50_l50();
    // No `sum_len == 0` guard: seqkit emits `NaN` for an all-zero-length
    // input (0.0/0.0), and guarding it would be defensive programming.
    let gc_percent = if matches!(seq_type, SeqType::Protein) {
        0.0
    } else {
        sum_gc as f64 * 100.0 / sum_len as f64
    };
    let avg_qual = if qual.err_sum == 0.0 {
        0.0
    } else {
        -10.0 * (qual.err_sum / sum_len as f64).log10()
    };
    ExtendedStats {
        q1: ls.q1(),
        q2: ls.q2(),
        q3: ls.q3(),
        sum_gap,
        n50,
        l50,
        q20_percent: round2(qual.q20 as f64 * 100.0 / sum_len as f64),
        q30_percent: round2(qual.q30 as f64 * 100.0 / sum_len as f64),
        avg_qual,
        gc_percent,
        sum_n,
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    fn write_fq(s: &str) -> tempfile::NamedTempFile {
        use std::io::Write;
        let mut f = tempfile::Builder::new()
            .suffix(".fq")
            .tempfile()
            .expect("tempfile");
        f.write_all(s.as_bytes()).expect("write");
        f
    }

    #[test]
    fn empty_input_errors() {
        let f = write_fq("");
        let err = compute_stats(f.path(), &Config::default()).unwrap_err();
        assert!(matches!(err, RsomicsError::InvalidInput(_)));
    }

    #[test]
    fn fasta_input_is_rejected() {
        let f = write_fq(">a\nACGT\n");
        let err = compute_stats(f.path(), &Config::default()).unwrap_err();
        assert!(matches!(err, RsomicsError::InvalidInput(_)));
    }

    // Locked against `seqkit stats -a -T` v2.13.0 on the identical input
    // (see tests/compat.rs for the live cross-check):
    // 3 reads, lens 10/6/12 → sum 28, avg 9.3; 26 of 28 bases q40, 2 q2.
    #[test]
    fn matches_seqkit_black_box_probe() {
        let f = write_fq(
            "@r1\nACGTACGTNN\n+\nIIIIIIII##\n@r2\nACGTAC\n+\nIIIIII\n@r3\nGGCCGGCCGGCC\n+\nIIIIIIIIIIII\n",
        );
        let cfg = Config {
            extended: true,
            ..Config::default()
        };
        let s = compute_stats(f.path(), &cfg).unwrap();
        assert_eq!(s.num_seqs, 3);
        assert_eq!(s.sum_len, 28);
        assert_eq!(s.min_len, 6);
        assert_eq!(s.max_len, 12);
        assert!((s.avg_len - 28.0 / 3.0).abs() < 1e-9);
        let e = s.extended.expect("--all");
        assert_eq!(e.n50, 10);
        assert_eq!(e.l50, 2);
        assert_eq!(e.sum_n, 2);
        // 26/28*100 = 92.857 → round2 → 92.86 → rendered %.0f = 93
        assert!((e.q20_percent - 92.86).abs() < 1e-9);
        assert!((e.q30_percent - 92.86).abs() < 1e-9);
        // errSum = 26·1e-4 + 2·10^-0.2 = 1.26452 ; -10·log10(1.26452/28)
        assert!((e.avg_qual - 13.45).abs() < 0.01);
        assert!((e.gc_percent - 19.0 * 100.0 / 28.0).abs() < 1e-9);
    }

    #[test]
    fn type_is_first_record_only() {
        let dna_then_protein =
            write_fq("@a\nACGTACGTACGT\n+\nIIIIIIIIIIII\n@b\nMEEPSILQRTWY\n+\nIIIIIIIIIIII\n");
        let s = compute_stats(dna_then_protein.path(), &Config::default()).unwrap();
        assert_eq!(s.seq_type, SeqType::Dna);

        let protein_then_dna =
            write_fq("@a\nMEEPSILQRTWY\n+\nIIIIIIIIIIII\n@b\nACGTACGTACGT\n+\nIIIIIIIIIIII\n");
        let s = compute_stats(protein_then_dna.path(), &Config::default()).unwrap();
        assert_eq!(s.seq_type, SeqType::Protein);
    }

    #[test]
    fn empty_first_record_is_unlimit() {
        let f = write_fq("@z\n\n+\n\n@b\nACGTACGT\n+\nIIIIIIII\n");
        let s = compute_stats(f.path(), &Config::default()).unwrap();
        assert_eq!(s.seq_type, SeqType::Unlimit);
    }

    #[test]
    fn wrong_encoding_offset_fails_loud() {
        let f = write_fq("@r\nACGT\n+\n!!!!\n");
        let cfg = Config {
            qual_offset: 64,
            ..Config::default()
        };
        let err = compute_stats(f.path(), &cfg).unwrap_err();
        assert!(matches!(err, RsomicsError::InvalidInput(_)));
    }
}
