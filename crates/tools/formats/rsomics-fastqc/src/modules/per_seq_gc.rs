// u64 read counts → f64 only at the final distribution/deviation stage.
#![allow(clippy::cast_precision_loss)]

use std::fmt::Write as _;

use super::{ModuleStatus, QcModule, Record};

/// `FastQC` "Per sequence GC content" — distribution of per-read GC% (0..=100)
/// against a theoretical normal centered on the *observed modal* GC. WARN if
/// the summed absolute deviation from that normal exceeds 15% of reads,
/// FAIL if it exceeds 30% (clean-room `FastQC` contract).
///
/// `FastQC` derives the reference curve from the observed data (modal GC as
/// the mean); the exact curve shape is calibrated against the `FastQC` binary
/// at the black-box compat step. The pass/warn/fail decision (the summed
/// deviation thresholds) is implemented per the documented contract here.
pub struct PerSeqGc {
    /// Observed read count per integer GC% bucket, 0..=100.
    obs: [u64; 101],
    total: u64,
}

impl PerSeqGc {
    #[must_use]
    pub fn new() -> Self {
        Self {
            obs: [0; 101],
            total: 0,
        }
    }

    /// Theoretical read count per GC% bucket: a normal centered on the
    /// observed mode with the observed standard deviation, scaled to the
    /// read total.
    fn theoretical(&self) -> [f64; 101] {
        let mut theo = [0.0_f64; 101];
        if self.total == 0 {
            return theo;
        }
        let total = self.total as f64;
        let mode = (0..=100).max_by_key(|&i| self.obs[i]).unwrap_or(50) as f64;
        let mut var = 0.0_f64;
        for (g, &c) in self.obs.iter().enumerate() {
            let d = g as f64 - mode;
            var += d * d * (c as f64);
        }
        let sigma = (var / total).sqrt().max(1.0);
        let norm = 1.0 / (sigma * (2.0 * std::f64::consts::PI).sqrt());
        for (g, t) in theo.iter_mut().enumerate() {
            let z = (g as f64 - mode) / sigma;
            *t = total * norm * (-0.5 * z * z).exp();
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
            dev += (self.obs[g] as f64 - t).abs();
        }
        // Each differing read is counted on both the observed and the
        // theoretical side, so halve to get the fraction of reads.
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
        let pct = (gc as f64 / rec.seq.len() as f64 * 100.0).round();
        // pct ∈ 0.0..=100.0 after round(); fits usize and is a valid index into obs[101]
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let pct_idx = pct as usize;
        self.obs[pct_idx] += 1;
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
        for (g, &c) in self.obs.iter().enumerate() {
            if c > 0 {
                let _ = writeln!(out, "{g}\t{c}");
            }
        }
    }
}
