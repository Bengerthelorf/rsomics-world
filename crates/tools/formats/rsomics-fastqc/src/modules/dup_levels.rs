// u64 counts → f64 only at the final percentage stage.
#![allow(clippy::cast_precision_loss)]

use std::collections::HashMap;
use std::fmt::Write as _;

use super::{ModuleStatus, QcModule, Record};

// FastQC also thins by sequence size beyond the read-index limit, so displayed curves may differ; status threshold is unaffected
pub struct DuplicationLevels {
    counts: HashMap<Vec<u8>, u64>,
    seen_reads: u64,
}

const TRACK_LIMIT: u64 = 100_000;
const KEY_TRUNC_OVER: usize = 75;
const KEY_LEN: usize = 50;

impl DuplicationLevels {
    #[must_use]
    pub fn new() -> Self {
        Self {
            counts: HashMap::new(),
            seen_reads: 0,
        }
    }

    fn key(seq: &[u8]) -> &[u8] {
        if seq.len() > KEY_TRUNC_OVER {
            &seq[..KEY_LEN]
        } else {
            seq
        }
    }

    fn distinct_total(&self) -> (u64, u64) {
        let distinct = self.counts.len() as u64;
        let total: u64 = self.counts.values().sum();
        (distinct, total)
    }

    fn non_unique_fraction(&self) -> f64 {
        let (distinct, total) = self.distinct_total();
        if total == 0 {
            0.0
        } else {
            1.0 - distinct as f64 / total as f64
        }
    }

    fn level_label(n: u64) -> &'static str {
        match n {
            0 => "0",
            1 => "1",
            2 => "2",
            3 => "3",
            4 => "4",
            5 => "5",
            6 => "6",
            7 => "7",
            8 => "8",
            9 => "9",
            10..=49 => ">10",
            50..=99 => ">50",
            100..=499 => ">100",
            500..=999 => ">500",
            1000..=4999 => ">1k",
            5000..=9999 => ">5k",
            _ => ">10k",
        }
    }
}

impl Default for DuplicationLevels {
    fn default() -> Self {
        Self::new()
    }
}

impl QcModule for DuplicationLevels {
    fn name(&self) -> &'static str {
        "Sequence Duplication Levels"
    }

    fn process(&mut self, rec: &Record) {
        self.seen_reads += 1;
        let key = Self::key(rec.seq);
        if let Some(c) = self.counts.get_mut(key) {
            *c += 1;
        } else if self.seen_reads <= TRACK_LIMIT {
            self.counts.insert(key.to_vec(), 1);
        }
    }

    fn finalize(&mut self) {}

    fn status(&self) -> ModuleStatus {
        let f = self.non_unique_fraction();
        if f > 0.50 {
            ModuleStatus::Fail
        } else if f > 0.20 {
            ModuleStatus::Warn
        } else {
            ModuleStatus::Pass
        }
    }

    fn write_data(&self, out: &mut String) {
        let (distinct, total) = self.distinct_total();
        let dedup_pct = if total == 0 {
            0.0
        } else {
            distinct as f64 / total as f64 * 100.0
        };
        let _ = writeln!(out, "#Total Deduplicated Percentage\t{dedup_pct:.6}");
        out.push_str("#Duplication Level\tPercentage of deduplicated\tPercentage of total\n");

        let mut by_level: HashMap<&'static str, (u64, u64)> = HashMap::new();
        for &c in self.counts.values() {
            let e = by_level.entry(Self::level_label(c)).or_insert((0, 0));
            e.0 += 1;
            e.1 += c;
        }
        for label in [
            "1", "2", "3", "4", "5", "6", "7", "8", "9", ">10", ">50", ">100", ">500", ">1k",
            ">5k", ">10k",
        ] {
            let (d, reads) = by_level.get(label).copied().unwrap_or((0, 0));
            let pd = if distinct == 0 {
                0.0
            } else {
                d as f64 / distinct as f64 * 100.0
            };
            let pt = if total == 0 {
                0.0
            } else {
                reads as f64 / total as f64 * 100.0
            };
            let _ = writeln!(out, "{label}\t{pd:.6}\t{pt:.6}");
        }
    }
}
