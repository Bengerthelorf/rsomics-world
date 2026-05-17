// fastp arithmetic uses int; values bounded by read length (≲hundreds of bp) — casts cannot overflow in practice.
#![allow(
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss
)]

// offset < 0 when insert is shorter than the read (adapter read-through).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OverlapResult {
    pub overlapped: bool,
    pub offset: isize,
    pub overlap_len: usize,
    pub diff: usize,
}

const NOT_OVERLAPPED: OverlapResult = OverlapResult {
    overlapped: false,
    offset: 0,
    overlap_len: 0,
    diff: 0,
};

// fastp Sequence complement: non-ACGTU (incl. N, lowercase) maps to N.
#[must_use]
pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}

#[must_use]
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

// fastp OverlapAnalysis::analyze, ported from AfterQC. diff_percent is a fraction (fastp --overlap_diff_percent_limit 20 → 0.2).
#[must_use]
pub fn analyze(
    r1: &[u8],
    r2: &[u8],
    diff_limit: usize,
    overlap_require: usize,
    diff_percent: f64,
) -> OverlapResult {
    let rcr2 = revcomp(r2);
    let len1 = r1.len() as isize;
    let len2 = rcr2.len() as isize;
    let str2 = &rcr2;
    let ovr = overlap_require as isize;
    let complete_compare_require: usize = 50;

    // Accept: fastp's `diff <= L || (diff > L && i > CCR)`; `diff > L` dropped (implied by `!(diff <= L)`).
    // Forward: rc(r2) slid right; when len1 <= ovr this loop is vacuously skipped and the reverse pass alone runs.
    let mut offset: isize = 0;
    while offset < len1 - ovr {
        let overlap_len = (len1 - offset).min(len2) as usize;
        let eff_limit = diff_limit.min((overlap_len as f64 * diff_percent) as usize);
        let mut diff = 0usize;
        let mut i = 0usize;
        while i < overlap_len {
            if r1[(offset as usize) + i] != str2[i] {
                diff += 1;
                if diff > eff_limit && i < complete_compare_require {
                    break;
                }
            }
            i += 1;
        }
        if diff <= eff_limit || i > complete_compare_require {
            return OverlapResult {
                overlapped: true,
                offset,
                overlap_len,
                diff,
            };
        }
        offset += 1;
    }

    // negative offset: insert is shorter than the read (adapter read-through).
    offset = 0;
    while offset > -(len2 - ovr) {
        let overlap_len = len1.min(len2 - offset.abs()) as usize;
        let eff_limit = diff_limit.min((overlap_len as f64 * diff_percent) as usize);
        let mut diff = 0usize;
        let mut i = 0usize;
        while i < overlap_len {
            if r1[i] != str2[(-offset as usize) + i] {
                diff += 1;
                if diff > eff_limit && i < complete_compare_require {
                    break;
                }
            }
            i += 1;
        }
        if diff <= eff_limit || i > complete_compare_require {
            return OverlapResult {
                overlapped: true,
                offset,
                overlap_len,
                diff,
            };
        }
        offset -= 1;
    }

    NOT_OVERLAPPED
}

// len1/len2: segment lengths for the fastp merged_<len1>_<len2> read name.
pub struct Merged {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub len1: usize,
    pub len2: usize,
}

// Caller must construct the name: fastp uses `r1.name + " merged_<len1>_<len2>"`.
#[must_use]
pub fn merge(
    r1_seq: &[u8],
    r1_qual: &[u8],
    r2_seq: &[u8],
    r2_qual: &[u8],
    ov: &OverlapResult,
) -> Option<Merged> {
    if !ov.overlapped {
        return None;
    }
    let ol = ov.overlap_len;
    let len1 = ol + ov.offset.max(0) as usize;
    let len2 = if ov.offset > 0 { r2_seq.len() - ol } else { 0 };

    let mut seq = r1_seq[..len1].to_vec();
    let mut qual = r1_qual[..len1].to_vec();
    if ov.offset > 0 {
        // fastp appends revcomp(r2)[ol..ol+len2]; with len2==r2.len()-ol that segment equals revcomp(r2[..len2]), so we extend from r2 directly — no full-length rc buffer (fastp allocates one per pair).
        seq.extend(r2_seq[..len2].iter().rev().map(|&b| complement(b)));
        qual.extend(r2_qual[..len2].iter().rev().copied());
    }
    Some(Merged {
        seq,
        qual,
        len1,
        len2,
    })
}

#[must_use]
pub const fn num2qual(n: u8) -> u8 {
    n + 33
}

// fastp BaseCorrector::correctByOverlapAnalysis: ≥Q30 vs ≤Q14 mismatch → rewrite to high-quality base.
pub fn correct(
    r1_seq: &mut [u8],
    r1_qual: &mut [u8],
    r2_seq: &mut [u8],
    r2_qual: &mut [u8],
    ov: &OverlapResult,
) -> usize {
    if ov.diff == 0 || !ov.overlapped {
        return 0;
    }
    let ol = ov.overlap_len;
    let start1 = ov.offset.max(0) as usize;
    let start2 = r2_seq.len() - (-ov.offset).max(0) as usize - 1;
    let good = num2qual(30);
    let bad = num2qual(14);

    let mut corrected = 0usize;
    for i in 0..ol {
        let p1 = start1 + i;
        let p2 = start2 - i;
        if r1_seq[p1] != complement(r2_seq[p2]) {
            if r1_qual[p1] >= good && r2_qual[p2] <= bad {
                r2_seq[p2] = complement(r1_seq[p1]);
                r2_qual[p2] = r1_qual[p1];
                corrected += 1;
            } else if r2_qual[p2] >= good && r1_qual[p1] <= bad {
                r1_seq[p1] = complement(r2_seq[p2]);
                r1_qual[p1] = r2_qual[p2];
                corrected += 1;
            }
        }
    }
    corrected
}

#[cfg(test)]
mod tests {
    use super::*;

    // Verbatim from fastp v0.20.1 src/overlapanalysis.cpp OverlapAnalysis::test().
    #[test]
    fn analyze_matches_fastp_test_vector() {
        let r1 = b"CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGC";
        let r2 = b"ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGTCC";
        let ov = analyze(r1, r2, 2, 30, 0.2);
        assert!(ov.overlapped);
        assert_eq!(ov.offset, 10);
        assert_eq!(ov.overlap_len, 79);
        assert_eq!(ov.diff, 1);
    }

    // Verbatim from fastp v0.20.1 src/basecorrector.cpp BaseCorrector::test().
    #[test]
    fn correct_matches_fastp_test_vector() {
        let mut r1s = b"TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCACGGGG".to_vec();
        let mut r1q = b"EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EEEEE".to_vec();
        let mut r2s = b"AAAAAAAAAACCCCGGGGAAAATTTTAAAATTGGGGGGGGGGTGGGGGGGGGGGGG".to_vec();
        let mut r2q = b"EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEEE".to_vec();
        let ov = analyze(&r1s, &r2s, 5, 30, 0.2);
        correct(&mut r1s, &mut r1q, &mut r2s, &mut r2q, &ov);
        assert_eq!(
            r1s,
            b"TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG"
        );
        assert_eq!(
            r2s,
            b"AAAAAAAAAACCCCGGGGAAAATTTTAAAATTGGGGGGGGGGGGGGGGGGGGGGGG"
        );
        assert_eq!(
            r1q,
            b"EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
        );
        assert_eq!(
            r2q,
            b"EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
        );
    }

    #[test]
    fn no_overlap_returns_not_overlapped() {
        let ov = analyze(
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
            5,
            30,
            0.2,
        );
        assert!(!ov.overlapped);
    }
}
