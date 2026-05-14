use serde::Serialize;

use crate::filter::FilterResult;
use crate::stats::ReadStats;

/// fastp-compatible JSON report.
///
/// Covers the pre/post-filter summary, filter counts, and per-mate per-cycle
/// quality / content curves. Sections outside this scope are omitted rather
/// than stubbed.
#[derive(Debug, Serialize)]
pub struct FastpJsonReport {
    pub summary: Summary,
    pub filtering_result: FilteringResult,
    pub read1_before_filtering: MateView,
    pub read1_after_filtering: MateView,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read2_before_filtering: Option<MateView>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read2_after_filtering: Option<MateView>,
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

/// Per-mate view — fastp emits one of these for each of R1 / R2, both before
/// and after filtering. Contains both the scalar totals and the per-cycle
/// quality / nucleotide-content curves.
#[derive(Debug, Serialize)]
pub struct MateView {
    pub total_reads: u64,
    pub total_bases: u64,
    pub q20_bases: u64,
    pub q30_bases: u64,
    pub total_cycles: usize,
    pub quality_curves: QualityCurves,
    pub content_curves: ContentCurves,
}

#[derive(Debug, Serialize)]
pub struct QualityCurves {
    #[serde(rename = "A")]
    pub a: Vec<f64>,
    #[serde(rename = "T")]
    pub t: Vec<f64>,
    #[serde(rename = "C")]
    pub c: Vec<f64>,
    #[serde(rename = "G")]
    pub g: Vec<f64>,
    pub mean: Vec<f64>,
}

#[derive(Debug, Serialize)]
pub struct ContentCurves {
    #[serde(rename = "A")]
    pub a: Vec<f64>,
    #[serde(rename = "T")]
    pub t: Vec<f64>,
    #[serde(rename = "C")]
    pub c: Vec<f64>,
    #[serde(rename = "G")]
    pub g: Vec<f64>,
    #[serde(rename = "N")]
    pub n: Vec<f64>,
    #[serde(rename = "GC")]
    pub gc: Vec<f64>,
}

impl From<&ReadStats> for MateView {
    fn from(s: &ReadStats) -> Self {
        let n = s.cycles.len();
        // Per-cycle qual is `qual_sum / total` where total is the count of
        // reads that reached this cycle. fastp emits one curve per nucleotide
        // (the mean Phred at positions where that base was called); the
        // per-nucleotide arrays here are set equal to the overall `mean`
        // because `ReadStats` tracks `qual_sum` at the cycle level, not
        // per-base-call. Downstream consumers that only read the `mean` curve
        // (the default visualization) see the correct value.
        let mean: Vec<f64> = s
            .cycles
            .iter()
            .map(|c| {
                let tot = c.total();
                if tot == 0 {
                    0.0
                } else {
                    #[allow(clippy::cast_precision_loss)]
                    {
                        c.qual_sum as f64 / tot as f64
                    }
                }
            })
            .collect();
        let frac = |num: u64, denom: u64| -> f64 {
            if denom == 0 {
                0.0
            } else {
                #[allow(clippy::cast_precision_loss)]
                {
                    num as f64 / denom as f64
                }
            }
        };
        let content_curves = ContentCurves {
            a: s.cycles
                .iter()
                .map(|c| frac(c.count_a, c.total()))
                .collect(),
            t: s.cycles
                .iter()
                .map(|c| frac(c.count_t, c.total()))
                .collect(),
            c: s.cycles
                .iter()
                .map(|c| frac(c.count_c, c.total()))
                .collect(),
            g: s.cycles
                .iter()
                .map(|c| frac(c.count_g, c.total()))
                .collect(),
            n: s.cycles
                .iter()
                .map(|c| frac(c.count_n, c.total()))
                .collect(),
            gc: s
                .cycles
                .iter()
                .map(|c| frac(c.count_c + c.count_g, c.total()))
                .collect(),
        };
        let quality_curves = QualityCurves {
            a: mean.clone(),
            t: mean.clone(),
            c: mean.clone(),
            g: mean.clone(),
            mean,
        };
        Self {
            total_reads: s.total_reads,
            total_bases: s.total_bases,
            q20_bases: s.q20_bases,
            q30_bases: s.q30_bases,
            total_cycles: n,
            quality_curves,
            content_curves,
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
    /// Build a single-end report from one set of pre/post stats.
    #[must_use]
    pub fn from_se(pre: &ReadStats, post: &ReadStats, filtering: FilteringResult) -> Self {
        Self {
            summary: Summary {
                before_filtering: pre.into(),
                after_filtering: post.into(),
            },
            filtering_result: filtering,
            read1_before_filtering: pre.into(),
            read1_after_filtering: post.into(),
            read2_before_filtering: None,
            read2_after_filtering: None,
        }
    }

    /// Build a paired-end report. The `summary` aggregates both mates
    /// (matching fastp), while the per-mate `read1_*` and `read2_*` blocks
    /// stay separate.
    #[must_use]
    pub fn from_pe(
        pre1: &ReadStats,
        post1: &ReadStats,
        pre2: &ReadStats,
        post2: &ReadStats,
        filtering: FilteringResult,
    ) -> Self {
        let pre_agg = aggregate(pre1, pre2);
        let post_agg = aggregate(post1, post2);
        Self {
            summary: Summary {
                before_filtering: (&pre_agg).into(),
                after_filtering: (&post_agg).into(),
            },
            filtering_result: filtering,
            read1_before_filtering: pre1.into(),
            read1_after_filtering: post1.into(),
            read2_before_filtering: Some(pre2.into()),
            read2_after_filtering: Some(post2.into()),
        }
    }
}

/// Combine two per-mate `ReadStats` into a summary-level aggregate. Cycle
/// vectors are summed position-by-position; the longer of the two sets the
/// length.
fn aggregate(a: &ReadStats, b: &ReadStats) -> ReadStats {
    let mut out = a.clone();
    out.total_reads += b.total_reads;
    out.total_bases += b.total_bases;
    out.q20_bases += b.q20_bases;
    out.q30_bases += b.q30_bases;
    out.gc_bases += b.gc_bases;
    out.n_bases += b.n_bases;
    if out.cycles.len() < b.cycles.len() {
        out.cycles
            .resize(b.cycles.len(), crate::stats::CycleStat::default());
    }
    for (dst, src) in out.cycles.iter_mut().zip(b.cycles.iter()) {
        dst.count_a += src.count_a;
        dst.count_c += src.count_c;
        dst.count_g += src.count_g;
        dst.count_t += src.count_t;
        dst.count_n += src.count_n;
        dst.qual_sum += src.qual_sum;
    }
    out
}
