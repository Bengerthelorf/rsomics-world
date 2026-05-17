// u64 read counts → f64 only at the final distribution/deviation stage.
#![allow(clippy::cast_precision_loss)]

use std::fmt::Write as _;

use super::{ModuleStatus, QcModule, Record};

pub struct PerSeqGc {
    // FastQC splits each read's GC% linearly across adjacent integer buckets; rounding to one bin would alias into a comb and inflate deviation
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
            // mean ∈ 0.0..=100.0 ⇒ round() is a valid index into theo[101]
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
        // Σtheo == Σobs, so Σ|diff| double-counts; halving gives total-variation distance
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
        // pct ∈ 0.0..=100.0, so floor/ceil are valid indices into obs[101]
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
        // FastQC emits all 101 buckets with one decimal (counts are fractional from the linear split)
        for (g, &c) in self.obs.iter().enumerate() {
            let _ = writeln!(out, "{g}\t{c:.1}");
        }
    }
}
