use std::collections::BTreeMap;
use std::fmt::Write as _;

use super::{ModuleStatus, QcModule, Record};

pub struct SeqLengthDistribution {
    counts: BTreeMap<u64, u64>,
}

impl SeqLengthDistribution {
    #[must_use]
    pub fn new() -> Self {
        Self {
            counts: BTreeMap::new(),
        }
    }
}

impl Default for SeqLengthDistribution {
    fn default() -> Self {
        Self::new()
    }
}

impl QcModule for SeqLengthDistribution {
    fn name(&self) -> &'static str {
        "Sequence Length Distribution"
    }

    fn process(&mut self, rec: &Record) {
        *self.counts.entry(rec.seq.len() as u64).or_insert(0) += 1;
    }

    fn finalize(&mut self) {}

    fn status(&self) -> ModuleStatus {
        if self.counts.contains_key(&0) {
            ModuleStatus::Fail
        } else if self.counts.len() > 1 {
            ModuleStatus::Warn
        } else {
            ModuleStatus::Pass
        }
    }

    fn write_data(&self, out: &mut String) {
        out.push_str("#Length\tCount\n");
        for (len, count) in &self.counts {
            let _ = writeln!(out, "{len}\t{count}");
        }
    }
}
