use std::fmt::Write as _;

use super::{ModuleStatus, QcModule, Record};

// FastQC bins distant positions into groups for long reads; we emit per-position — status is exact, display layout differs for very long reads
pub struct PerBaseNContent {
    n: Vec<u64>,
    total: Vec<u64>,
}

impl PerBaseNContent {
    #[must_use]
    pub fn new() -> Self {
        Self {
            n: Vec::new(),
            total: Vec::new(),
        }
    }

    fn max_n_fraction(&self) -> f64 {
        let mut worst = 0.0_f64;
        for (i, &t) in self.total.iter().enumerate() {
            if t == 0 {
                continue;
            }
            #[allow(clippy::cast_precision_loss)]
            let frac = self.n[i] as f64 / t as f64;
            if frac > worst {
                worst = frac;
            }
        }
        worst
    }
}

impl Default for PerBaseNContent {
    fn default() -> Self {
        Self::new()
    }
}

impl QcModule for PerBaseNContent {
    fn name(&self) -> &'static str {
        "Per base N content"
    }

    fn process(&mut self, rec: &Record) {
        if rec.seq.len() > self.n.len() {
            self.n.resize(rec.seq.len(), 0);
            self.total.resize(rec.seq.len(), 0);
        }
        for (i, &b) in rec.seq.iter().enumerate() {
            self.total[i] += 1;
            if b == b'N' || b == b'n' {
                self.n[i] += 1;
            }
        }
    }

    fn finalize(&mut self) {}

    fn status(&self) -> ModuleStatus {
        let worst = self.max_n_fraction();
        if worst > 0.20 {
            ModuleStatus::Fail
        } else if worst > 0.05 {
            ModuleStatus::Warn
        } else {
            ModuleStatus::Pass
        }
    }

    fn write_data(&self, out: &mut String) {
        out.push_str("#Base\tN-Count\n");
        for (i, &t) in self.total.iter().enumerate() {
            #[allow(clippy::cast_precision_loss)]
            let pct = if t == 0 {
                0.0
            } else {
                self.n[i] as f64 * 100.0 / t as f64
            };
            let _ = writeln!(out, "{}\t{:.6}", i + 1, pct);
        }
    }
}
