use std::fmt::Write as _;

use super::{ModuleStatus, QcModule, Record};

/// `FastQC` "Per sequence quality scores" — distribution of each read's mean
/// quality. WARN if the most frequent mean quality (the mode) is below 27,
/// FAIL if below 20 (clean-room `FastQC` contract).
///
/// The per-read mean is taken over raw quality bytes (offset-independent);
/// the Phred offset is resolved in `finalize` from the lowest byte seen
/// (Phred+33 unless the lowest byte is ≥ 64), matching `FastQC`'s encoding
/// inference, then the histogram is shifted into Phred space.
pub struct PerSeqQuality {
    raw_mean_hist: [u64; 256],
    min_byte: u8,
    qhist: Vec<u64>,
    offset: u8,
}

impl PerSeqQuality {
    #[must_use]
    pub fn new() -> Self {
        Self {
            raw_mean_hist: [0; 256],
            min_byte: u8::MAX,
            qhist: Vec::new(),
            offset: 33,
        }
    }

    fn mode_quality(&self) -> Option<usize> {
        let mut best = None;
        let mut best_count = 0;
        for (q, &c) in self.qhist.iter().enumerate() {
            if c > best_count {
                best_count = c;
                best = Some(q);
            }
        }
        best
    }
}

impl Default for PerSeqQuality {
    fn default() -> Self {
        Self::new()
    }
}

impl QcModule for PerSeqQuality {
    fn name(&self) -> &'static str {
        "Per sequence quality scores"
    }

    fn process(&mut self, rec: &Record) {
        if rec.qual.is_empty() {
            return;
        }
        let mut sum: u64 = 0;
        for &q in rec.qual {
            sum += u64::from(q);
            self.min_byte = self.min_byte.min(q);
        }
        #[allow(clippy::cast_precision_loss)]
        let mean = (sum as f64 / rec.qual.len() as f64).round();
        // mean of bytes in 0..=255 ⇒ rounded mean ∈ 0.0..=255.0, fits usize
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let mean_idx = mean as usize;
        self.raw_mean_hist[mean_idx] += 1;
    }

    fn finalize(&mut self) {
        self.offset = if self.min_byte < 64 { 33 } else { 64 };
        let off = usize::from(self.offset);
        self.qhist = vec![0; 256 - off];
        for (raw, &count) in self.raw_mean_hist.iter().enumerate() {
            if count == 0 {
                continue;
            }
            let q = raw.saturating_sub(off);
            self.qhist[q] += count;
        }
    }

    fn status(&self) -> ModuleStatus {
        match self.mode_quality() {
            Some(m) if m < 20 => ModuleStatus::Fail,
            Some(m) if m < 27 => ModuleStatus::Warn,
            _ => ModuleStatus::Pass,
        }
    }

    fn write_data(&self, out: &mut String) {
        out.push_str("#Quality\tCount\n");
        for (q, &c) in self.qhist.iter().enumerate() {
            if c > 0 {
                let _ = writeln!(out, "{q}\t{c}");
            }
        }
    }
}
