// Hamming-distance only — fastp's `--allow_gap_overlap_trimming` path is not
// implemented. Trim geometry: negative offset = adapter present; non-negative
// = no adapter trim needed.

/// `diff_percent_limit` is clamped to `[0.0, 1.0]` at `sanitised` construction.
#[derive(Debug, Clone, Copy)]
pub struct OverlapConfig {
    pub overlap_require: usize,
    pub diff_limit: usize,
    pub diff_percent_limit: f32,
}

impl OverlapConfig {
    #[must_use]
    pub fn sanitised(overlap_require: usize, diff_limit: usize, diff_percent_limit: f32) -> Self {
        let pct = if diff_percent_limit.is_nan() {
            0.0
        } else {
            diff_percent_limit.clamp(0.0, 1.0)
        };
        Self {
            overlap_require,
            diff_limit,
            diff_percent_limit: pct,
        }
    }
}

impl Default for OverlapConfig {
    fn default() -> Self {
        Self {
            overlap_require: 30,
            diff_limit: 5,
            diff_percent_limit: 0.20,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct OverlapResult {
    pub overlapped: bool,
    /// Positive = R2-RC starts inside R1; negative = adapter-present geometry.
    pub offset: i64,
    pub overlap_len: usize,
    pub diff: usize,
}

impl OverlapResult {
    #[must_use]
    pub const fn none() -> Self {
        Self {
            overlapped: false,
            offset: 0,
            overlap_len: 0,
            diff: 0,
        }
    }
}

#[must_use]
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        out.push(complement(b));
    }
    out
}

#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        other => other,
    }
}

#[must_use]
pub fn analyze(r1: &[u8], r2_rc: &[u8], cfg: OverlapConfig) -> OverlapResult {
    let len1 = r1.len();
    let len2 = r2_rc.len();
    if len1 < cfg.overlap_require || len2 < cfg.overlap_require {
        return OverlapResult::none();
    }

    // Forward scan: offset in [0, len1 - overlap_require].
    let max_fwd = len1 - cfg.overlap_require;
    for offset in 0..=max_fwd {
        let overlap_len = (len1 - offset).min(len2);
        let budget = budget_for(overlap_len, cfg);
        let diff = count_mismatches_bounded(
            &r1[offset..offset + overlap_len],
            &r2_rc[..overlap_len],
            budget,
        );
        if diff <= budget {
            return OverlapResult {
                overlapped: true,
                offset: i64::try_from(offset).unwrap_or(i64::MAX),
                overlap_len,
                diff,
            };
        }
    }

    // Reverse scan: offset in [-1, -(len2 - overlap_require)].
    let max_rev = len2 - cfg.overlap_require;
    for shift in 1..=max_rev {
        let overlap_len = len1.min(len2 - shift);
        let budget = budget_for(overlap_len, cfg);
        let diff = count_mismatches_bounded(
            &r1[..overlap_len],
            &r2_rc[shift..shift + overlap_len],
            budget,
        );
        if diff <= budget {
            return OverlapResult {
                overlapped: true,
                offset: -i64::try_from(shift).unwrap_or(i64::MAX),
                overlap_len,
                diff,
            };
        }
    }

    OverlapResult::none()
}

#[inline]
fn budget_for(overlap_len: usize, cfg: OverlapConfig) -> usize {
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    let by_pct = (overlap_len as f32 * cfg.diff_percent_limit) as usize;
    cfg.diff_limit.min(by_pct)
}

#[inline]
fn count_mismatches_bounded(a: &[u8], b: &[u8], limit: usize) -> usize {
    debug_assert_eq!(a.len(), b.len());
    let mut diff = 0usize;
    for (x, y) in a.iter().zip(b.iter()) {
        if !x.eq_ignore_ascii_case(y) {
            diff += 1;
            if diff > limit {
                return diff;
            }
        }
    }
    diff
}

/// `front_trimmed1`/`front_trimmed2` are already-applied fixed-front trims —
/// required to keep overlap-derived lengths consistent with the original read frame.
#[must_use]
pub fn trim_lengths(
    ov: OverlapResult,
    len1: usize,
    len2: usize,
    front_trimmed1: usize,
    front_trimmed2: usize,
) -> Option<(usize, usize)> {
    if !ov.overlapped || ov.offset >= 0 {
        return None;
    }
    let new_len1 = len1.min(ov.overlap_len + front_trimmed2);
    let new_len2 = len2.min(ov.overlap_len + front_trimmed1);
    Some((new_len1, new_len2))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rc_basic() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"ACGTN"), b"NACGT");
    }

    #[test]
    fn no_overlap_when_reads_too_short() {
        let r1 = b"ACGT";
        let r2_rc = b"ACGT";
        let ov = analyze(r1, r2_rc, OverlapConfig::default());
        assert!(!ov.overlapped);
    }

    #[test]
    fn perfect_full_overlap_at_offset_zero() {
        // 30 bp R1, R2-RC equals R1 → full overlap at offset 0, no adapter.
        let r1 = b"ACGTACGTACGTACGTACGTACGTACGTAC";
        let r2_rc = r1;
        let ov = analyze(r1, r2_rc, OverlapConfig::default());
        assert!(ov.overlapped);
        assert_eq!(ov.offset, 0);
        assert_eq!(ov.overlap_len, 30);
    }

    #[test]
    fn adapter_present_negative_offset_case() {
        // Insert (20 bp) shared between R1 and R2-RC. Adapter (5 bp) past
        // each end. Geometry:
        //   R1     = insert ++ adapter1                       (25 bp)
        //   R2     = rc(insert) ++ adapter2                   (25 bp)
        //   R2-RC  = rc(adapter2) ++ insert                   (25 bp)
        // Reverse scan with shift=5 makes r1[..20] == r2_rc[5..25].
        let insert: &[u8] = b"AAACCCCGGGGTTTTCCAAT";
        let adapter1 = b"GAATC";
        let adapter2 = b"GGTTA";
        let mut r1 = insert.to_vec();
        r1.extend_from_slice(adapter1);
        let mut r2 = reverse_complement(insert);
        r2.extend_from_slice(adapter2);
        let r2_rc = reverse_complement(&r2);
        let cfg = OverlapConfig {
            overlap_require: 10,
            ..OverlapConfig::default()
        };
        let ov = analyze(&r1, &r2_rc, cfg);
        assert!(ov.overlapped, "{ov:?}");
        assert_eq!(ov.offset, -5, "{ov:?}");
        assert_eq!(ov.overlap_len, 20);
    }

    #[test]
    fn trim_lengths_skips_positive_offset() {
        let ov = OverlapResult {
            overlapped: true,
            offset: 5,
            overlap_len: 25,
            diff: 0,
        };
        assert_eq!(trim_lengths(ov, 30, 30, 0, 0), None);
    }

    #[test]
    fn trim_lengths_negative_offset_returns_truncations() {
        let ov = OverlapResult {
            overlapped: true,
            offset: -10,
            overlap_len: 20,
            diff: 0,
        };
        assert_eq!(trim_lengths(ov, 30, 30, 0, 0), Some((20, 20)));
    }
}
