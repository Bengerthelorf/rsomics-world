use serde::Serialize;

use crate::filter::FilterResult;
use crate::stats::ReadStats;

/// fastp-compatible JSON report subset.
///
/// fastp's full schema is ~40 sections; this covers the four most-used: the
/// pre/post-filter summaries and the filtering-result counts. Additional
/// sections (`adapter_cutting`, `duplication`, `insert_size`, etc.) will land
/// as the corresponding pipeline stages do.
#[derive(Debug, Serialize)]
pub struct FastpJsonReport {
    pub summary: Summary,
    pub filtering_result: FilteringResult,
}

#[derive(Debug, Serialize)]
pub struct Summary {
    pub before_filtering: ReadStatsView,
    pub after_filtering: ReadStatsView,
}

#[derive(Debug, Serialize)]
pub struct ReadStatsView {
    pub total_reads: u64,
    pub total_bases: u64,
    pub q20_bases: u64,
    pub q30_bases: u64,
    pub q20_rate: f64,
    pub q30_rate: f64,
    pub gc_content: f64,
}

impl From<&ReadStats> for ReadStatsView {
    fn from(s: &ReadStats) -> Self {
        Self {
            total_reads: s.total_reads,
            total_bases: s.total_bases,
            q20_bases: s.q20_bases,
            q30_bases: s.q30_bases,
            q20_rate: s.q20_rate(),
            q30_rate: s.q30_rate(),
            gc_content: s.gc_content(),
        }
    }
}

#[derive(Debug, Default, Serialize)]
pub struct FilteringResult {
    pub passed_filter_reads: u64,
    pub low_quality_reads: u64,
    #[serde(rename = "too_many_N_reads")]
    pub too_many_n_reads: u64,
    pub too_short_reads: u64,
}

impl FilteringResult {
    /// Increment the appropriate counter for one read's filter outcome.
    pub fn record(&mut self, outcome: FilterResult) {
        match outcome {
            FilterResult::Pass => self.passed_filter_reads += 1,
            FilterResult::LowQuality => self.low_quality_reads += 1,
            FilterResult::TooManyN => self.too_many_n_reads += 1,
            FilterResult::TooShort => self.too_short_reads += 1,
        }
    }
}

impl FastpJsonReport {
    #[must_use]
    pub fn from_stats(pre: &ReadStats, post: &ReadStats, filtering: FilteringResult) -> Self {
        Self {
            summary: Summary {
                before_filtering: pre.into(),
                after_filtering: post.into(),
            },
            filtering_result: filtering,
        }
    }
}
