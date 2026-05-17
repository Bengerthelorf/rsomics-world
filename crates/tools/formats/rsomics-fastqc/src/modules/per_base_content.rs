// u64 base counts → f64 percentages only at the final ratio stage; counts
// fit the 52-bit mantissa for any real FASTQ.
#![allow(clippy::cast_precision_loss)]

use std::fmt::Write as _;

use super::{ModuleStatus, QcModule, Record};

fn pct(n: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        n as f64 / total as f64 * 100.0
    }
}

/// `FastQC` "Per base sequence content" — %A/C/G/T at each position. WARN if
/// |%A−%T| or |%G−%C| exceeds 10 at any position; FAIL if it exceeds 20
/// (clean-room `FastQC` contract). `FastQC` emits columns in G,A,T,C order.
pub struct PerBaseContent {
    a: Vec<u64>,
    c: Vec<u64>,
    g: Vec<u64>,
    t: Vec<u64>,
    total: Vec<u64>,
}

impl PerBaseContent {
    #[must_use]
    pub fn new() -> Self {
        Self {
            a: Vec::new(),
            c: Vec::new(),
            g: Vec::new(),
            t: Vec::new(),
            total: Vec::new(),
        }
    }

    /// Largest (|A−T|, |G−C|) percentage gap over all positions.
    fn worst_gap(&self) -> f64 {
        let mut worst = 0.0_f64;
        for i in 0..self.total.len() {
            let tot = self.total[i];
            if tot == 0 {
                continue;
            }
            let at = (pct(self.a[i], tot) - pct(self.t[i], tot)).abs();
            let gc = (pct(self.g[i], tot) - pct(self.c[i], tot)).abs();
            worst = worst.max(at).max(gc);
        }
        worst
    }
}

impl Default for PerBaseContent {
    fn default() -> Self {
        Self::new()
    }
}

impl QcModule for PerBaseContent {
    fn name(&self) -> &'static str {
        "Per base sequence content"
    }

    fn process(&mut self, rec: &Record) {
        if rec.seq.len() > self.total.len() {
            self.a.resize(rec.seq.len(), 0);
            self.c.resize(rec.seq.len(), 0);
            self.g.resize(rec.seq.len(), 0);
            self.t.resize(rec.seq.len(), 0);
            self.total.resize(rec.seq.len(), 0);
        }
        for (i, &b) in rec.seq.iter().enumerate() {
            self.total[i] += 1;
            match b {
                b'A' | b'a' => self.a[i] += 1,
                b'C' | b'c' => self.c[i] += 1,
                b'G' | b'g' => self.g[i] += 1,
                b'T' | b't' => self.t[i] += 1,
                _ => {}
            }
        }
    }

    fn finalize(&mut self) {}

    fn status(&self) -> ModuleStatus {
        let w = self.worst_gap();
        if w > 20.0 {
            ModuleStatus::Fail
        } else if w > 10.0 {
            ModuleStatus::Warn
        } else {
            ModuleStatus::Pass
        }
    }

    fn write_data(&self, out: &mut String) {
        out.push_str("#Base\tG\tA\tT\tC\n");
        for i in 0..self.total.len() {
            let tot = self.total[i];
            if tot == 0 {
                continue;
            }
            let _ = writeln!(
                out,
                "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                i + 1,
                pct(self.g[i], tot),
                pct(self.a[i], tot),
                pct(self.t[i], tot),
                pct(self.c[i], tot),
            );
        }
    }
}
