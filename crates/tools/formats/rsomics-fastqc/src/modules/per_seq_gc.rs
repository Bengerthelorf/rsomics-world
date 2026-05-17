// u64 read counts → f64 only at the final distribution/deviation stage.
#![allow(clippy::cast_precision_loss)]

use std::fmt::Write as _;

use super::{ModuleStatus, QcModule, Record};

/// `FastQC` "Per sequence GC content" — the per-read GC% distribution
/// (0..=100) compared to a theoretical normal fitted to the observed data.
/// WARN if the deviation from that normal exceeds 15% of reads, FAIL if it
/// exceeds 30% (clean-room `FastQC` contract).
///
/// The reference is a normal with the observed mean and standard deviation,
/// evaluated on the discrete 0..=100 GC% support and renormalised so it
/// sums to the read total. Comparing like-with-like (both distributions
/// over the same discrete support, same total) is what makes a clean
/// unimodal library deviate ≈0 (PASS) while a bimodal/contaminated library
/// deviates sharply. The deviation is the total-variation distance
/// `Σ|obs−theo| / 2 / total` ∈ [0,1] — the fraction of reads that would
/// have to move to turn the observed distribution into the reference.
pub struct PerSeqGc {
    /// Observed read weight per integer GC% bucket, 0..=100. Each read's
    /// real-valued GC% is split linearly between its two adjacent integer
    /// buckets (`FastQC`'s method), so the distribution is smooth — bucketing
    /// to a single rounded bin would alias into a comb for read lengths not
    /// a multiple of 100 and spuriously inflate the deviation.
    obs: [f64; 101],
    total: u64,
}

impl PerSeqGc {
    #[must_use]
    pub fn new() -> Self {
        Self {
            obs: [0.0; 101],
            total: 0,
        }
    }

    /// Reference read count per GC% bucket: a normal with the observed mean
    /// and SD over the discrete 0..=100 support, renormalised so the total
    /// equals the read count. A zero-variance library (every read the same
    /// GC%) collapses to a single spike at that bucket.
    fn theoretical(&self) -> [f64; 101] {
        let mut theo = [0.0_f64; 101];
        if self.total == 0 {
            return theo;
        }
        let total = self.total as f64;
        let mut mean = 0.0_f64;
        for (g, &c) in self.obs.iter().enumerate() {
            mean += g as f64 * c;
        }
        mean /= total;
        let mut var = 0.0_f64;
        for (g, &c) in self.obs.iter().enumerate() {
            let d = g as f64 - mean;
            var += d * d * c;
        }
        let sigma = (var / total).sqrt();
        if sigma == 0.0 {
            // Degenerate: all reads share one GC% — the reference is a
            // single spike there, matching the observed spike exactly.
            // mean ∈ 0.0..=100.0 ⇒ round() is a valid index into theo[101].
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let spike = mean.round() as usize;
            theo[spike] = total;
            return theo;
        }
        let mut raw = [0.0_f64; 101];
        let mut sum_raw = 0.0_f64;
        for (g, r) in raw.iter_mut().enumerate() {
            let z = (g as f64 - mean) / sigma;
            *r = (-0.5 * z * z).exp();
            sum_raw += *r;
        }
        for (g, t) in theo.iter_mut().enumerate() {
            *t = raw[g] * total / sum_raw;
        }
        theo
    }

    fn deviation_fraction(&self) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        let theo = self.theoretical();
        let mut dev = 0.0_f64;
        for (g, &t) in theo.iter().enumerate() {
            dev += (self.obs[g] - t).abs();
        }
        // Σtheo == Σobs == total, so the absolute difference double-counts
        // every moved read (once as a deficit, once as a surplus); halving
        // gives the total-variation distance = fraction of reads moved.
        dev / 2.0 / self.total as f64
    }
}

impl Default for PerSeqGc {
    fn default() -> Self {
        Self::new()
    }
}

impl QcModule for PerSeqGc {
    fn name(&self) -> &'static str {
        "Per sequence GC content"
    }

    fn process(&mut self, rec: &Record) {
        if rec.seq.is_empty() {
            return;
        }
        let gc = rec
            .seq
            .iter()
            .filter(|&&b| matches!(b, b'G' | b'C' | b'g' | b'c'))
            .count();
        // FastQC splits each read's real-valued GC% linearly across the two
        // adjacent integer buckets (yielding a smooth distribution); a
        // single rounded bin would alias into a comb. pct ∈ 0.0..=100.0 so
        // floor/ceil are valid indices into obs[101].
        let pct = gc as f64 / rec.seq.len() as f64 * 100.0;
        let lo = pct.floor();
        let hi = pct.ceil();
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let (lo_i, hi_i) = (lo as usize, hi as usize);
        if lo_i == hi_i {
            self.obs[lo_i] += 1.0;
        } else {
            self.obs[lo_i] += hi - pct;
            self.obs[hi_i] += pct - lo;
        }
        self.total += 1;
    }

    fn finalize(&mut self) {}

    fn status(&self) -> ModuleStatus {
        let f = self.deviation_fraction();
        if f > 0.30 {
            ModuleStatus::Fail
        } else if f > 0.15 {
            ModuleStatus::Warn
        } else {
            ModuleStatus::Pass
        }
    }

    fn write_data(&self, out: &mut String) {
        out.push_str("#GC Content\tCount\n");
        // FastQC emits every bucket 0..=100 with one decimal (counts are
        // fractional because each read is split across two buckets).
        for (g, &c) in self.obs.iter().enumerate() {
            let _ = writeln!(out, "{g}\t{c:.1}");
        }
    }
}
