// one shared implementation for every rsomics-*-stats tool — re-deriving the math per format would silently break --all --tabular byte-agreement with seqkit stats -a -T
// lengths/counts fit f64's 52-bit mantissa for any real genome
#![allow(clippy::cast_precision_loss)]

use serde::Serialize;

#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq)]
pub enum SeqType {
    #[serde(rename = "DNA")]
    Dna,
    #[serde(rename = "RNA")]
    Rna,
    #[serde(rename = "Protein")]
    Protein,
    // seqkit's catch-all is "Unlimit" not "Other" — byte-equality needs this exact string
    #[serde(rename = "Unlimit")]
    Unlimit,
}

impl SeqType {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Dna => "DNA",
            Self::Rna => "RNA",
            Self::Protein => "Protein",
            Self::Unlimit => "Unlimit",
        }
    }
}

// seqkit AlphabetGuessSeqLengthThreshold; first record only — accumulating across records diverges
pub const DEFAULT_ALPHABET_GUESS_LEN: usize = 10_000;

// dedup needles: overlapping classes (e.g. b"GCgc") must not double-count
#[must_use]
pub fn count_any_of(haystack: &[u8], needles: &[u8]) -> u64 {
    let mut seen = [false; 256];
    let mut total: u64 = 0;
    for &n in needles {
        if seen[n as usize] {
            continue;
        }
        seen[n as usize] = true;
        // per-needle bytecount SIMD (4 passes) beats a scalar single-pass LUT — NEON/AVX2 recoups the bandwidth (host-M: 83 ms vs 53 ms, chr22 fixture)
        total += bytecount::count(haystack, n) as u64;
    }
    total
}

// seqkit seq.GuessAlphabet — pass the first record's prefix only
#[must_use]
pub fn classify(sample: &[u8]) -> SeqType {
    if sample.is_empty() {
        return SeqType::Unlimit;
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
            _ => return SeqType::Unlimit,
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

// port of shenwei356/bio length-stats.go — L50 = unique-length-bucket count (not records), matches seqkit --tabular --all
pub struct LengthStats {
    counts: Vec<(u64, u64)>,
    sum: u64,
    count: u64,
}

impl LengthStats {
    #[must_use]
    pub fn new(mut lengths: Vec<u64>) -> Self {
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
        // callers pass i_med_l < self.count, so the last bucket's acc always terminates the loop
        unreachable!(
            "LengthStats::get_value: i_med_l={i_med_l} not bracketed in counts (count={}, flag={flag}, prev={prev})",
            self.count
        )
    }

    #[must_use]
    pub fn q2(&self) -> f64 {
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

    #[must_use]
    pub fn q1(&self) -> f64 {
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

    #[must_use]
    pub fn q3(&self) -> f64 {
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

    #[must_use]
    pub fn n50_l50(&self) -> (u64, u64) {
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

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

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
    fn n50_three_contigs() {
        let ls = LengthStats::new(vec![4, 6, 8]);
        assert_eq!(ls.n50_l50(), (6, 2));
    }

    #[test]
    fn single_length_bucket() {
        let ls = LengthStats::new(vec![7, 7, 7]);
        assert_eq!(ls.q1(), 7.0);
        assert_eq!(ls.q2(), 7.0);
        assert_eq!(ls.q3(), 7.0);
        assert_eq!(ls.n50_l50(), (7, 1));
    }

    #[test]
    fn empty_is_zero() {
        let ls = LengthStats::new(vec![]);
        assert_eq!(ls.q2(), 0.0);
        assert_eq!(ls.n50_l50(), (0, 0));
    }

    #[test]
    fn count_any_of_dedupes_overlapping_classes() {
        assert_eq!(count_any_of(b"GCgcAT", b"GCgc"), 4);
        assert_eq!(count_any_of(b"NNnnX", b"Nn"), 4);
    }

    #[test]
    fn classify_rna_protein_dna_unlimit() {
        assert_eq!(classify(b"ACGU"), SeqType::Rna);
        assert_eq!(classify(b"MEEPSILQRT"), SeqType::Protein);
        assert_eq!(classify(b"ACGTN"), SeqType::Dna);
        assert_eq!(classify(b""), SeqType::Unlimit);
        assert_eq!(classify(b"@#$"), SeqType::Unlimit);
    }
}
