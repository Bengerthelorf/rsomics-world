//! Poly-X tail trimming. Two scan modes:
//!
//! - [`find_polyx_3p`] — forced-base scan (default poly-G).
//! - [`find_dominant_polyx_3p`] — A/C/G/T simultaneous count, dominant
//!   base post-scan.
//!
//! Ports of `PolyX::trimPolyG` and `PolyX::trimPolyX` from fastp's
//! `polyx.cpp`. The rate-based mismatch budget absorbs interspersed
//! non-target bases inside a long run — an isolated G at the 5'-most
//! edge of a poly-G stretch still shifts the trim point left.

use std::num::NonZeroUsize;

#[derive(Debug, Clone, Copy)]
pub struct PolyXConfig {
    pub base: u8,
    pub min_len: usize,
    pub max_mismatches: usize,
    pub mismatch_per_bases: NonZeroUsize,
}

impl Default for PolyXConfig {
    fn default() -> Self {
        Self {
            base: b'G',
            min_len: 10,
            max_mismatches: 5,
            mismatch_per_bases: NonZeroUsize::new(8).expect("8 is nonzero"),
        }
    }
}

impl PolyXConfig {
    #[must_use]
    pub fn poly_g() -> Self {
        Self::default()
    }

    #[must_use]
    pub fn for_base(base: u8) -> Self {
        Self {
            base,
            ..Self::default()
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct DominantPolyXResult {
    /// Dominant base of the trimmed tail (uppercase ASCII).
    pub base: u8,
    /// 0-based offset where the trim happens — keep `seq[..trim_at]`.
    pub trim_at: usize,
}

/// Forced-base scan. Returns the 0-based trim offset, or `None` if no
/// qualifying tail of length `≥ cfg.min_len` is found.
#[must_use]
pub fn find_polyx_3p(seq: &[u8], cfg: PolyXConfig) -> Option<usize> {
    let rlen = seq.len();
    if rlen == 0 {
        return None;
    }
    let target = cfg.base.to_ascii_uppercase();
    let mut mismatch: usize = 0;
    let mut last_target_pos: Option<usize> = None;
    let mut i: usize = 0;
    let divisor = cfg.mismatch_per_bases.get();
    while i < rlen {
        let b = seq[rlen - i - 1].to_ascii_uppercase();
        if b == target {
            last_target_pos = Some(rlen - i - 1);
        } else {
            mismatch += 1;
        }
        let allowed = (i + 1) / divisor;
        if mismatch > cfg.max_mismatches || (mismatch > allowed && i + 1 >= cfg.min_len) {
            break;
        }
        i += 1;
    }
    if i >= cfg.min_len {
        last_target_pos
    } else {
        None
    }
}

/// Dominant-base scan. Counts A/C/G/T simultaneously walking from the 3'
/// end; on stop, picks the most-represented base as the polyX target and
/// returns the trim offset at its last-occurrence position. `N` increments
/// all four counters (fastp parity).
#[must_use]
pub fn find_dominant_polyx_3p(seq: &[u8], cfg: PolyXConfig) -> Option<DominantPolyXResult> {
    let rlen = seq.len();
    if rlen == 0 {
        return None;
    }
    // counts[0]=A, [1]=C, [2]=G, [3]=T
    let mut counts: [usize; 4] = [0; 4];
    let mut i: usize = 0;
    let divisor = cfg.mismatch_per_bases.get();
    while i < rlen {
        let b = seq[rlen - i - 1].to_ascii_uppercase();
        match b {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            b'N' => {
                for c in &mut counts {
                    *c += 1;
                }
            }
            _ => {}
        }
        let cmp = i + 1;
        let allowed = cfg.max_mismatches.min(cmp / divisor);
        let mut any_within_budget = false;
        for &c in &counts {
            if cmp - c <= allowed {
                any_within_budget = true;
                break;
            }
        }
        if !any_within_budget && (i >= divisor || cmp + 1 >= cfg.min_len) {
            break;
        }
        i += 1;
    }
    if i + 1 < cfg.min_len {
        return None;
    }
    let mut dom_idx = 0usize;
    for j in 1..4 {
        if counts[j] > counts[dom_idx] {
            dom_idx = j;
        }
    }
    let dom_base = b"ACGT"[dom_idx];
    // `i` counts scanned bytes (0..=rlen); back-search uses
    // `seq[rlen - pos - 1]`, so `pos` must stay within `[0, rlen-1]`.
    // The natural-exit path (loop ran to `i == rlen`) needs the clamp.
    let mut pos = i.min(rlen - 1);
    while pos > 0 && seq[rlen - pos - 1].to_ascii_uppercase() != dom_base {
        pos -= 1;
    }
    if seq[rlen - pos - 1].to_ascii_uppercase() != dom_base {
        return None;
    }
    Some(DominantPolyXResult {
        base: dom_base,
        trim_at: rlen - pos - 1,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_polyg_means_no_trim() {
        let seq = b"ACGTACGTACGTACGTACGTACGT";
        assert_eq!(find_polyx_3p(seq, PolyXConfig::default()), None);
    }

    #[test]
    fn polyg_at_3prime_trims_at_run_start() {
        let seq = b"ACGTACGTACGTACGTACGTGGGGGGGGGGGG";
        assert_eq!(find_polyx_3p(seq, PolyXConfig::default()), Some(18));
    }

    #[test]
    fn polyg_below_min_len_is_not_trimmed() {
        let seq = b"ACGTACGTACGTACGTACGTGGGGG";
        assert_eq!(find_polyx_3p(seq, PolyXConfig::default()), None);
    }

    #[test]
    fn mismatches_within_budget_are_tolerated() {
        let seq = b"ACGTACGTACGTGGGGGAGGGG";
        let cfg = PolyXConfig {
            min_len: 9,
            max_mismatches: 1,
            mismatch_per_bases: NonZeroUsize::new(8).unwrap(),
            ..Default::default()
        };
        assert_eq!(find_polyx_3p(seq, cfg), Some(12));
    }

    #[test]
    fn lowercase_g_counts() {
        let seq = b"ACGTACGTACGTACGTACGTgggggggggggg";
        assert_eq!(find_polyx_3p(seq, PolyXConfig::default()), Some(18));
    }

    #[test]
    fn polya_with_for_base() {
        let seq = b"ACGTACGTACGTACGTACGTAAAAAAAAAAAA";
        let cfg = PolyXConfig::for_base(b'A');
        assert_eq!(find_polyx_3p(seq, cfg), Some(20));
    }

    #[test]
    fn dominant_polyx_detects_polya() {
        let seq = b"GCTAGCTAGCTAGCTAGCTAAAAAAAAAAAA";
        let r = find_dominant_polyx_3p(seq, PolyXConfig::default()).unwrap();
        assert_eq!(r.base, b'A');
    }

    #[test]
    fn dominant_polyx_detects_polyc() {
        let seq = b"GATGCATGCATGCATGCATGCCCCCCCCCCCC";
        let r = find_dominant_polyx_3p(seq, PolyXConfig::default()).unwrap();
        assert_eq!(r.base, b'C');
    }

    #[test]
    fn dominant_polyx_no_dominant_tail_returns_none() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        assert!(find_dominant_polyx_3p(seq, PolyXConfig::default()).is_none());
    }

    #[test]
    fn dominant_polyx_iupac_ambiguity_no_panic() {
        // IUPAC ambiguity codes (R/Y/S/W/K/M/B/D/H/V) are not in the
        // A/C/G/T/N count arms — counts stay zero. With a loose budget the
        // scan reaches `i == rlen`; the back-search must still be in-bounds.
        let seq = b"RYSWKMBDHVRYSWKMBDHVRYSWKMBDHVRY";
        let cfg = PolyXConfig {
            base: b'G',
            min_len: seq.len(),
            max_mismatches: usize::MAX,
            mismatch_per_bases: NonZeroUsize::new(1).unwrap(),
        };
        assert!(find_dominant_polyx_3p(seq, cfg).is_none());
    }
}
