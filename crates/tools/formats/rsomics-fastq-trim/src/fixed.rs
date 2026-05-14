//! Fixed-length 5' / 3' trimming. Drops a fixed prefix and/or suffix
//! regardless of content. Matches fastp's `--trim_front{1,2}` and
//! `--trim_tail{1,2}` flags.

/// Fixed-trim configuration. Per-mate so PE callers can use different
/// front/tail trims on R1 vs R2 — matches fastp's separate `_1` / `_2`
/// flag pairs.
#[derive(Debug, Clone, Copy, Default)]
pub struct FixedTrimConfig {
    pub trim_front: usize,
    pub trim_tail: usize,
}

/// Compute the half-open `[start, end)` slice of the record after applying
/// the fixed-length cut. Returns `None` if `front + tail >= seq.len()`,
/// i.e. the cut would leave a zero-length read — the caller is expected
/// to discard such records.
#[must_use]
pub fn apply_fixed(seq_len: usize, cfg: FixedTrimConfig) -> Option<(usize, usize)> {
    let total = cfg.trim_front.checked_add(cfg.trim_tail)?;
    if total >= seq_len {
        return None;
    }
    Some((cfg.trim_front, seq_len - cfg.trim_tail))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_trim_returns_full_range() {
        assert_eq!(apply_fixed(20, FixedTrimConfig::default()), Some((0, 20)));
    }

    #[test]
    fn front_only() {
        let cfg = FixedTrimConfig {
            trim_front: 6,
            trim_tail: 0,
        };
        assert_eq!(apply_fixed(20, cfg), Some((6, 20)));
    }

    #[test]
    fn tail_only() {
        let cfg = FixedTrimConfig {
            trim_front: 0,
            trim_tail: 4,
        };
        assert_eq!(apply_fixed(20, cfg), Some((0, 16)));
    }

    #[test]
    fn both_combine() {
        let cfg = FixedTrimConfig {
            trim_front: 6,
            trim_tail: 4,
        };
        assert_eq!(apply_fixed(20, cfg), Some((6, 16)));
    }

    #[test]
    fn over_trim_returns_none() {
        let cfg = FixedTrimConfig {
            trim_front: 10,
            trim_tail: 10,
        };
        assert_eq!(apply_fixed(20, cfg), None);
    }
}
