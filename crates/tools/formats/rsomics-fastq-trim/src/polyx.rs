//! Poly-X tail trimming. Generalises fastp's poly-G handling (the default
//! for `NextSeq` / `NovaSeq` 2-color chemistry) to any single-base run.
//!
//! Verbatim port of `PolyX::trimPolyG` from fastp's `polyx.cpp`, with
//! the target base parameterised so the same scan handles poly-A / -T /
//! -C / -G alike. The algorithm:
//!
//! - Walk from the 3' end one base at a time.
//! - Track the last position where the target base was seen (`first_g_pos`
//!   in fastp's naming — for us "last target seen").
//! - Allow mismatches up to two simultaneous caps:
//!   - a hard cap (fastp default 5)
//!   - a rate cap of `floor(scanned / divisor)` (fastp uses divisor 8 =
//!     "one mismatch per 8 bases").
//! - Stop when either cap is exceeded *and* we've scanned ≥ `min_len`
//!   bases; trim at the last-seen-target position.

/// Poly-X configuration.
#[derive(Debug, Clone, Copy)]
pub struct PolyXConfig {
    pub base: u8,
    pub min_len: usize,
    pub max_mismatches: usize,
    pub mismatch_per_bases: usize,
}

impl Default for PolyXConfig {
    fn default() -> Self {
        Self {
            base: b'G',
            min_len: 10,
            max_mismatches: 5,
            mismatch_per_bases: 8,
        }
    }
}

impl PolyXConfig {
    #[must_use]
    pub fn poly_g() -> Self {
        Self::default()
    }

    #[must_use]
    pub fn poly_x(base: u8) -> Self {
        Self {
            base,
            ..Self::default()
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct PolyXResult {
    pub trimmed: usize,
}

/// Return the 0-based offset at which the 3' poly-X tail starts, or
/// `None` if no qualifying tail is found. Trim point is the last
/// occurrence of `cfg.base` within the scanned tail — non-base
/// interspersed bases that fall within the allowed mismatch budget
/// stay inside the trimmed-off region (matches fastp's behaviour
/// where an isolated G at the 5'-most end of the poly-G run still
/// shifts the trim point left).
#[must_use]
pub fn find_polyx_3p(seq: &[u8], cfg: PolyXConfig) -> Option<usize> {
    let rlen = seq.len();
    if rlen == 0 {
        return None;
    }
    let target = cfg.base.to_ascii_uppercase();
    let mut mismatch: usize = 0;
    let mut last_target_pos: usize = rlen - 1;
    let mut i: usize = 0;
    let mut found_any_target = false;
    while i < rlen {
        let b = seq[rlen - i - 1].to_ascii_uppercase();
        if b == target {
            last_target_pos = rlen - i - 1;
            found_any_target = true;
        } else {
            mismatch += 1;
        }
        let allowed = (i + 1) / cfg.mismatch_per_bases;
        if mismatch > cfg.max_mismatches || (mismatch > allowed && i + 1 >= cfg.min_len) {
            break;
        }
        i += 1;
    }
    if i >= cfg.min_len && found_any_target {
        Some(last_target_pos)
    } else {
        None
    }
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
        // fastp's rate-based algorithm absorbs the isolated G at pos 18
        // into the run because the mismatch rate (1/14 bases) stays within
        // the 1-per-8 budget. Trim point is the last G observed.
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
            mismatch_per_bases: 8,
            ..Default::default()
        };
        assert_eq!(find_polyx_3p(seq, cfg), Some(12));
    }

    #[test]
    fn lowercase_g_counts() {
        let seq = b"ACGTACGTACGTACGTACGTgggggggggggg";
        // Same rate-based absorption as the uppercase test — isolated G
        // at pos 18 is included in the run.
        assert_eq!(find_polyx_3p(seq, PolyXConfig::default()), Some(18));
    }

    #[test]
    fn polya_with_base_param() {
        let seq = b"ACGTACGTACGTACGTACGTAAAAAAAAAAAA";
        let cfg = PolyXConfig::poly_x(b'A');
        assert_eq!(find_polyx_3p(seq, cfg), Some(20));
    }
}
