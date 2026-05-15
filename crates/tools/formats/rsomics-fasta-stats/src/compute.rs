// u64→f64 cast: lengths/counts fit in 52-bit mantissa for any real genome;
// cast only at final quartile/N50/percentage stage.
#![allow(clippy::cast_precision_loss)]

use std::path::Path;

use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

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

#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq)]
pub enum SeqType {
    #[serde(rename = "DNA")]
    Dna,
    #[serde(rename = "RNA")]
    Rna,
    #[serde(rename = "Protein")]
    Protein,
    #[serde(rename = "Other")]
    Other,
}

impl SeqType {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Dna => "DNA",
            Self::Rna => "RNA",
            Self::Protein => "Protein",
            Self::Other => "Other",
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

    // bytecount SIMD (4-pass) beats scalar single-pass: NEON/AVX2 recoups
    // bandwidth cost. (M2: scalar LUT 83 ms vs 53 ms on chr22 fixture.)
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

fn count_any_of(haystack: &[u8], needles: &[u8]) -> u64 {
    let mut seen = [false; 256];
    let mut total: u64 = 0;
    for &n in needles {
        if seen[n as usize] {
            continue;
        }
        seen[n as usize] = true;
        total += bytecount::count(haystack, n) as u64;
    }
    total
}

// Port of bio/util/length-stats.go. seqkit's L50 counts unique-length buckets,
// not records — reproduced so `--tabular --all` agrees with seqkit.
struct LengthStats {
    counts: Vec<(u64, u64)>, // (length, cumulative_count) sorted ascending, deduped
    sum: u64,
    count: u64,
}

impl LengthStats {
    fn new(mut lengths: Vec<u64>) -> Self {
        let sum: u64 = lengths.iter().sum();
        let count = lengths.len() as u64;
        lengths.sort_unstable();
        let mut counts: Vec<(u64, u64)> = Vec::new();
        let mut acc: u64 = 0;
        let mut i = 0;
        while i < lengths.len() {
            let v = lengths[i];
            let mut j = i;
            while j < lengths.len() && lengths[j] == v {
                j += 1;
            }
            acc += (j - i) as u64;
            counts.push((v, acc));
            i = j;
        }
        Self { counts, sum, count }
    }

    fn get_value(&self, even: bool, i_med_l: u64, i_med_r: u64) -> f64 {
        let mut flag = false;
        let mut prev: u64 = 0;
        for &(len, acc) in &self.counts {
            if flag {
                return (len + prev) as f64 / 2.0;
            }
            if acc > i_med_l {
                if even {
                    if acc > i_med_r {
                        return len as f64;
                    }
                    flag = true;
                    prev = len;
                } else {
                    return len as f64;
                }
            }
        }
        // Callers pass i_med_l < self.count, so the last bucket's acc always
        // terminates the loop. Reaching here means a caller broke the invariant.
        unreachable!(
            "LengthStats::get_value: i_med_l={i_med_l} not bracketed in counts (count={}, flag={flag}, prev={prev})",
            self.count
        )
    }

    fn q2(&self) -> f64 {
        if self.counts.is_empty() {
            return 0.0;
        }
        if self.counts.len() == 1 {
            return self.counts[0].0 as f64;
        }
        let even = self.count & 1 == 0;
        if even {
            let l = self.count / 2 - 1;
            let r = self.count / 2;
            self.get_value(true, l, r)
        } else {
            self.get_value(false, self.count / 2, 0)
        }
    }

    fn q1(&self) -> f64 {
        if self.counts.is_empty() {
            return 0.0;
        }
        if self.counts.len() == 1 {
            return self.counts[0].0 as f64;
        }
        let parent_even = self.count & 1 == 0;
        let n = if parent_even {
            self.count / 2
        } else {
            self.count.div_ceil(2)
        };
        let even = n & 1 == 0;
        if even {
            self.get_value(true, n / 2 - 1, n / 2)
        } else {
            self.get_value(false, n / 2, 0)
        }
    }

    fn q3(&self) -> f64 {
        if self.counts.is_empty() {
            return 0.0;
        }
        if self.counts.len() == 1 {
            return self.counts[0].0 as f64;
        }
        let parent_even = self.count & 1 == 0;
        let (n, mean) = if parent_even {
            (self.count / 2, self.count / 2)
        } else {
            (self.count.div_ceil(2), self.count / 2)
        };
        let even = n & 1 == 0;
        if even {
            self.get_value(true, n / 2 - 1 + mean, n / 2 + mean)
        } else {
            self.get_value(false, n / 2 + mean, 0)
        }
    }

    fn n50_l50(&self) -> (u64, u64) {
        if self.counts.is_empty() {
            return (0, 0);
        }
        if self.counts.len() == 1 {
            return (self.counts[0].0, 1);
        }
        let half = self.sum as f64 / 2.0;
        let mut sum_len: f64 = 0.0;
        let n = self.counts.len();
        for i in (0..n).rev() {
            let (len, acc) = self.counts[i];
            let prev_acc = if i == 0 { 0 } else { self.counts[i - 1].1 };
            let per_len_count = acc - prev_acc;
            sum_len += (len * per_len_count) as f64;
            if sum_len >= half {
                return (len, (n - i) as u64);
            }
        }
        (0, 0)
    }
}

fn classify(sample: &[u8]) -> SeqType {
    if sample.is_empty() {
        return SeqType::Other;
    }
    let mut has_t = false;
    let mut has_u = false;
    let mut has_protein_only = false;
    for &b in sample {
        let c = b.to_ascii_uppercase();
        match c {
            b'T' => has_t = true,
            b'U' => has_u = true,
            b'E' | b'F' | b'I' | b'L' | b'P' | b'Q' | b'Z' | b'X' | b'*' => {
                has_protein_only = true;
            }
            b'A' | b'C' | b'G' | b'N' | b'-' | b'.' | b' ' | b'\n' | b'\r' | b'R' | b'Y' | b'S'
            | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' => {}
            _ => return SeqType::Other,
        }
    }
    if has_protein_only {
        SeqType::Protein
    } else if has_u && !has_t {
        SeqType::Rna
    } else {
        SeqType::Dna
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
    fn quartiles_three_contigs_match_seqkit() {
        let ls = LengthStats::new(vec![4, 6, 8]);
        assert_eq!(ls.q1(), 5.0);
        assert_eq!(ls.q2(), 6.0);
        assert_eq!(ls.q3(), 7.0);
    }

    #[test]
    fn quartiles_one_to_nine_match_seqkit() {
        let ls = LengthStats::new((1u64..=9).collect());
        assert_eq!(ls.q1(), 3.0);
        assert_eq!(ls.q2(), 5.0);
        assert_eq!(ls.q3(), 7.0);
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
