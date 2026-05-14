/// Threshold configuration for read-level filtering.
///
/// Defaults match fastp's defaults so that pipelines using `rsomics-fastp` with
/// no overrides reject the same reads fastp would reject.
#[derive(Debug, Clone, Copy)]
pub struct FilterConfig {
    pub qualified_quality_phred: u8,
    pub unqualified_percent_limit: u32,
    pub length_required: usize,
    pub n_base_limit: usize,
}

impl Default for FilterConfig {
    fn default() -> Self {
        Self {
            qualified_quality_phred: 15,
            unqualified_percent_limit: 40,
            length_required: 15,
            n_base_limit: 5,
        }
    }
}

/// Per-read filter outcome. Ordered to match fastp's evaluation precedence:
/// length first (cheapest), then N count, then quality.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterResult {
    Pass,
    TooShort,
    TooManyN,
    LowQuality,
}

/// Classify a single read against `cfg`. Quality bytes are Sanger-encoded
/// (Phred + 33).
#[must_use]
pub fn classify(seq: &[u8], qual: &[u8], cfg: FilterConfig) -> FilterResult {
    if seq.len() < cfg.length_required {
        return FilterResult::TooShort;
    }

    let n_count = seq.iter().filter(|&&b| b == b'N' || b == b'n').count();
    if n_count > cfg.n_base_limit {
        return FilterResult::TooManyN;
    }

    let threshold = cfg.qualified_quality_phred + 33;
    let unqualified = qual.iter().filter(|&&q| q < threshold).count();
    let pct = (unqualified * 100) / qual.len().max(1);
    if pct > cfg.unqualified_percent_limit as usize {
        return FilterResult::LowQuality;
    }

    FilterResult::Pass
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pass_high_q() {
        let seq = b"ACGTACGTACGTACGT";
        let qual = b"IIIIIIIIIIIIIIII";
        assert_eq!(
            classify(seq, qual, FilterConfig::default()),
            FilterResult::Pass
        );
    }

    #[test]
    fn too_short_rejects_below_required_length() {
        let seq = b"ACGT";
        let qual = b"IIII";
        assert_eq!(
            classify(seq, qual, FilterConfig::default()),
            FilterResult::TooShort
        );
    }

    #[test]
    fn too_many_n_rejects_above_limit() {
        let seq = b"NNNNNNACGTACGTACGTACGT";
        let qual = b"IIIIIIIIIIIIIIIIIIIIII";
        assert_eq!(
            classify(seq, qual, FilterConfig::default()),
            FilterResult::TooManyN
        );
    }

    #[test]
    fn low_quality_rejects_when_unqualified_fraction_exceeds_limit() {
        // All bases at Q2 (`#`), 100% unqualified > 40%.
        let seq = b"ACGTACGTACGTACGT";
        let qual = b"################";
        assert_eq!(
            classify(seq, qual, FilterConfig::default()),
            FilterResult::LowQuality
        );
    }

    #[test]
    fn precedence_length_before_n() {
        // 4 N's in 4-bp read: too short fires first.
        let seq = b"NNNN";
        let qual = b"IIII";
        assert_eq!(
            classify(seq, qual, FilterConfig::default()),
            FilterResult::TooShort
        );
    }
}
