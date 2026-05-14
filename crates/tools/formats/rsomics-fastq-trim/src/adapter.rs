//! 3' adapter trimming via prefix-match scan.
//!
//! Scan each offset in the read, compare the suffix against the adapter
//! prefix within a mismatch budget, trim at the earliest qualifying match.
//! Covers SE and the R1-only PE fallback path; the PE overlap-detection
//! variant lives in [`crate::overlap`].
//!
//! Algorithm: Hamming-distance only — fastp's one-insertion and
//! one-deletion fallback phases are not implemented.

/// Adapter-trim configuration. Defaults match fastp's: `TruSeq` adapter
/// (R1 sense-strand prefix), 5 bp minimum compared length, 20% mismatch
/// budget.
#[derive(Debug, Clone)]
pub struct AdapterConfig {
    pub sequence: Vec<u8>,
    pub min_match_len: usize,
    pub max_mismatch_rate: f32,
}

impl AdapterConfig {
    #[must_use]
    pub fn illumina_truseq_r1() -> Self {
        Self {
            sequence: b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec(),
            min_match_len: 5,
            max_mismatch_rate: 0.2,
        }
    }

    #[must_use]
    pub fn illumina_truseq_r2() -> Self {
        Self {
            sequence: b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_vec(),
            min_match_len: 5,
            max_mismatch_rate: 0.2,
        }
    }
}

/// Trim outcome for one record.
#[derive(Debug, Clone, Copy, Default)]
pub struct AdapterResult {
    /// Number of bases trimmed off the 3' end (0 if no adapter found).
    pub trimmed: usize,
}

/// Return the 0-based offset where the read should be trimmed to remove
/// the adapter, or `None` if no signature found at the 3' end. The
/// earliest qualifying match wins.
///
/// "Qualifying match" = at least `min_match_len` bases compared and the
/// mismatch fraction across the compared region ≤ `max_mismatch_rate`.
#[must_use]
pub fn find_adapter_3p(seq: &[u8], cfg: &AdapterConfig) -> Option<usize> {
    let adapter = &cfg.sequence;
    if adapter.is_empty() || seq.len() < cfg.min_match_len {
        return None;
    }

    let max_start = seq.len().saturating_sub(cfg.min_match_len);
    for start in 0..=max_start {
        let cmp_len = (seq.len() - start).min(adapter.len());
        if cmp_len < cfg.min_match_len {
            continue;
        }
        let mismatches = seq[start..start + cmp_len]
            .iter()
            .zip(&adapter[..cmp_len])
            .filter(|(a, b)| !a.eq_ignore_ascii_case(b))
            .count();
        #[allow(clippy::cast_precision_loss)]
        let rate = mismatches as f32 / cmp_len as f32;
        if rate <= cfg.max_mismatch_rate {
            return Some(start);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_adapter_means_no_trim() {
        let seq = b"ACGTACGTACGTACGTACGTACGT";
        let cfg = AdapterConfig::illumina_truseq_r1();
        assert_eq!(find_adapter_3p(seq, &cfg), None);
    }

    #[test]
    fn perfect_adapter_at_3prime_is_trimmed() {
        let seq = b"ACGTACGTACGTACGTACGTAGATCGGAAGAG";
        let cfg = AdapterConfig::illumina_truseq_r1();
        assert_eq!(find_adapter_3p(seq, &cfg), Some(20));
    }

    #[test]
    fn partial_adapter_within_mismatch_budget() {
        let seq = b"ACGTACGTACGTACGTACGTAAATCGGAAGAG";
        let cfg = AdapterConfig::illumina_truseq_r1();
        assert_eq!(find_adapter_3p(seq, &cfg), Some(20));
    }

    #[test]
    fn too_few_bases_to_match_returns_none() {
        let seq = b"AGAT";
        let cfg = AdapterConfig::illumina_truseq_r1();
        assert_eq!(find_adapter_3p(seq, &cfg), None);
    }

    #[test]
    fn case_insensitive_matching() {
        let seq = b"ACGTACGTACGTACGTACGTagatcggaagag";
        let cfg = AdapterConfig::illumina_truseq_r1();
        assert_eq!(find_adapter_3p(seq, &cfg), Some(20));
    }
}
