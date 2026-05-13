//! `PolyG` tail trimming for 2-color-chemistry instruments (`NextSeq` / `NovaSeq`).
//!
//! When the cluster signal goes dark at the 3' end, the basecaller emits G
//! (G = no signal in 2-color chemistry) and assigns spuriously high quality
//! scores, so vanilla quality trimming misses it. We detect a run of G's
//! anchored at the 3' end within a small mismatch budget and trim it off
//! before the adapter / filter pipeline sees the read.

/// PolyG-trimming configuration. Defaults match fastp's defaults: 10 bp
/// minimum poly-G run.
#[derive(Debug, Clone, Copy)]
pub struct PolyGConfig {
    pub min_len: usize,
    pub max_mismatches: usize,
}

impl Default for PolyGConfig {
    fn default() -> Self {
        Self {
            min_len: 10,
            max_mismatches: 0,
        }
    }
}

/// Return the 0-based offset to trim the 3' poly-G tail at, or `None` if no
/// qualifying tail is found. "Qualifying" = at least `min_len` bases at the
/// 3' end where the count of non-G bases is ≤ `max_mismatches`.
///
/// The trim point is the earliest (5'-most) position from which all but at
/// most `max_mismatches` bases are G — i.e. we trim as aggressively as the
/// budget allows.
#[must_use]
pub fn find_polyg_3p(seq: &[u8], cfg: PolyGConfig) -> Option<usize> {
    if seq.len() < cfg.min_len {
        return None;
    }

    let max_start = seq.len() - cfg.min_len;
    let mut best: Option<usize> = None;
    for start in (0..=max_start).rev() {
        let tail = &seq[start..];
        let mismatches = tail.iter().filter(|&&b| !is_g(b)).count();
        if mismatches > cfg.max_mismatches {
            // Extending further left can only add more non-G bases.
            break;
        }
        best = Some(start);
    }
    best
}

fn is_g(b: u8) -> bool {
    b == b'G' || b == b'g'
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_polyg_means_no_trim() {
        let seq = b"ACGTACGTACGTACGTACGTACGT";
        assert_eq!(find_polyg_3p(seq, PolyGConfig::default()), None);
    }

    #[test]
    fn polyg_at_3prime_trims_at_run_start() {
        // 20 bp insert + 12 G tail.
        let seq = b"ACGTACGTACGTACGTACGTGGGGGGGGGGGG";
        assert_eq!(find_polyg_3p(seq, PolyGConfig::default()), Some(20));
    }

    #[test]
    fn polyg_below_min_len_is_not_trimmed() {
        // 5 G's < min_len 10.
        let seq = b"ACGTACGTACGTACGTACGTGGGGG";
        assert_eq!(find_polyg_3p(seq, PolyGConfig::default()), None);
    }

    #[test]
    fn mismatches_within_budget_are_tolerated() {
        // 10 G's with one A in the middle, budget = 1.
        let seq = b"ACGTACGTACGTGGGGGAGGGG";
        let cfg = PolyGConfig {
            min_len: 9,
            max_mismatches: 1,
        };
        // Greedy: the longest qualifying tail starts at the leftmost position
        // where the running non-G count is still ≤ 1.
        assert_eq!(find_polyg_3p(seq, cfg), Some(12));
    }

    #[test]
    fn lowercase_g_counts() {
        let seq = b"ACGTACGTACGTACGTACGTgggggggggggg";
        assert_eq!(find_polyg_3p(seq, PolyGConfig::default()), Some(20));
    }
}
