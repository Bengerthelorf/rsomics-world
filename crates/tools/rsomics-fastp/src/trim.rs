//! Adapter trimming.
//!
//! For single-end reads we use sequence-based 3' adapter matching: scan each
//! offset in the read, compare against the adapter prefix within a mismatch
//! budget, and trim at the earliest qualifying match. Paired-end overlap-based
//! adapter detection (fastp's higher-quality default) is a follow-up.

/// Adapter-trimming configuration. Defaults match fastp's defaults: `TruSeq`
/// adapter, 0.2 mismatch rate, 5 bp minimum match.
#[derive(Debug, Clone)]
pub struct AdapterConfig {
    pub sequence: Vec<u8>,
    pub min_match_len: usize,
    pub max_mismatch_rate: f32,
}

impl AdapterConfig {
    /// Illumina ``TruSeq`` adapter (sense-strand prefix, R1).
    #[must_use]
    pub fn illumina_truseq_r1() -> Self {
        Self {
            sequence: b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec(),
            min_match_len: 5,
            max_mismatch_rate: 0.2,
        }
    }
}

/// Return the 0-based offset where the read should be trimmed to remove the
/// adapter, or `None` if no adapter signature is found at the 3' end. The
/// earliest qualifying match wins so we trim as aggressively as the budget
/// allows.
///
/// "Qualifying match" = at least `min_match_len` bases compared (i.e. the
/// suffix has to be long enough), and the mismatch fraction across the
/// compared region is ≤ `max_mismatch_rate`.
#[must_use]
pub fn find_adapter_3p(seq: &[u8], cfg: &AdapterConfig) -> Option<usize> {
    let adapter = &cfg.sequence;
    if adapter.is_empty() || seq.len() < cfg.min_match_len {
        return None;
    }

    // Scan from earliest position (most aggressive trim) to latest.
    let max_start = seq.len().saturating_sub(cfg.min_match_len);
    for start in 0..=max_start {
        let cmp_len = (seq.len() - start).min(adapter.len());
        if cmp_len < cfg.min_match_len {
            continue;
        }
        let mismatches = seq[start..start + cmp_len]
            .iter()
            .zip(&adapter[..cmp_len])
            .filter(|(a, b)| a.eq_ignore_ascii_case(b).not())
            .count();
        // cast: cmp_len fits in u32 trivially (read length).
        #[allow(clippy::cast_precision_loss)]
        let rate = mismatches as f32 / cmp_len as f32;
        if rate <= cfg.max_mismatch_rate {
            return Some(start);
        }
    }
    None
}

// Local extension so the closure above reads cleanly.
trait NotExt {
    fn not(self) -> bool;
}
impl NotExt for bool {
    fn not(self) -> bool {
        !self
    }
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
        // 20 bp insert + adapter prefix.
        let seq = b"ACGTACGTACGTACGTACGTAGATCGGAAGAG";
        let cfg = AdapterConfig::illumina_truseq_r1();
        assert_eq!(find_adapter_3p(seq, &cfg), Some(20));
    }

    #[test]
    fn partial_adapter_within_mismatch_budget() {
        // 20 bp insert + adapter prefix with 1 mismatch in the first 5 bp.
        let seq = b"ACGTACGTACGTACGTACGTAAATCGGAAGAG";
        let cfg = AdapterConfig::illumina_truseq_r1();
        // 1/12 = 8% mismatch, under 20% budget.
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
