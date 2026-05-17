use std::fmt::Write as _;

use super::{ModuleStatus, QcModule, Record};

// FastQC bins distant positions into groups for long reads; we emit per-position — status is exact, display layout differs for very long reads
pub struct PerBaseQuality {
    hist: Vec<[u64; 256]>,
    min_byte: u8,
    offset: u8,
}

struct PosStats {
    mean: f64,
    median: f64,
    lower_q: f64,
    upper_q: f64,
    p10: f64,
    p90: f64,
}

impl PerBaseQuality {
    #[must_use]
    pub fn new() -> Self {
        Self {
            hist: Vec::new(),
            min_byte: u8::MAX,
            offset: 33,
        }
    }

    fn percentile(counts: &[u64; 256], offset: u8, total: u64, pct: f64) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let target = pct / 100.0 * total as f64;
        let off = usize::from(offset);
        let mut cum: u64 = 0;
        for (byte, &c) in counts.iter().enumerate() {
            if c == 0 {
                continue;
            }
            cum += c;
            #[allow(clippy::cast_precision_loss)]
            if cum as f64 >= target {
                return byte.saturating_sub(off) as f64;
            }
        }
        0.0
    }

    fn pos_stats(&self, pos: usize) -> Option<PosStats> {
        let counts = &self.hist[pos];
        let total: u64 = counts.iter().sum();
        if total == 0 {
            return None;
        }
        let off = usize::from(self.offset);
        let mut wsum: u64 = 0;
        for (byte, &c) in counts.iter().enumerate() {
            if c == 0 {
                continue;
            }
            wsum += byte.saturating_sub(off) as u64 * c;
        }
        #[allow(clippy::cast_precision_loss)]
        let mean = wsum as f64 / total as f64;
        Some(PosStats {
            mean,
            median: Self::percentile(counts, self.offset, total, 50.0),
            lower_q: Self::percentile(counts, self.offset, total, 25.0),
            upper_q: Self::percentile(counts, self.offset, total, 75.0),
            p10: Self::percentile(counts, self.offset, total, 10.0),
            p90: Self::percentile(counts, self.offset, total, 90.0),
        })
    }
}

impl Default for PerBaseQuality {
    fn default() -> Self {
        Self::new()
    }
}

impl QcModule for PerBaseQuality {
    fn name(&self) -> &'static str {
        "Per base sequence quality"
    }

    fn process(&mut self, rec: &Record) {
        if rec.qual.len() > self.hist.len() {
            self.hist.resize(rec.qual.len(), [0; 256]);
        }
        for (i, &q) in rec.qual.iter().enumerate() {
            self.hist[i][q as usize] += 1;
            self.min_byte = self.min_byte.min(q);
        }
    }

    fn finalize(&mut self) {
        self.offset = if self.min_byte < 64 { 33 } else { 64 };
    }

    fn status(&self) -> ModuleStatus {
        let mut worst = ModuleStatus::Pass;
        for pos in 0..self.hist.len() {
            let Some(s) = self.pos_stats(pos) else {
                continue;
            };
            if s.lower_q < 5.0 || s.median < 20.0 {
                return ModuleStatus::Fail;
            }
            if s.lower_q < 10.0 || s.median < 25.0 {
                worst = ModuleStatus::Warn;
            }
        }
        worst
    }

    fn write_data(&self, out: &mut String) {
        out.push_str(
            "#Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile\n",
        );
        for pos in 0..self.hist.len() {
            let Some(s) = self.pos_stats(pos) else {
                continue;
            };
            let _ = writeln!(
                out,
                "{}\t{:.6}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{:.1}",
                pos + 1,
                s.mean,
                s.median,
                s.lower_q,
                s.upper_q,
                s.p10,
                s.p90
            );
        }
    }
}
