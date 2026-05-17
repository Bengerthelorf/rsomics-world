// u64 counts → f64 only at the final percentage stage.
#![allow(clippy::cast_precision_loss)]

use std::fmt::Write as _;

use super::{ModuleStatus, QcModule, Record};

/// `FastQC` "Adapter Content" — cumulative % of reads in which each adapter
/// has appeared at or before each position (once seen in a read it stays
/// positive to the read's end). WARN if any adapter exceeds 5% at any
/// position, FAIL if it exceeds 10% (clean-room `FastQC` contract).
///
/// Adapter k-mers are `FastQC`'s public `adapter_list.txt` (12 bp each).
const ADAPTERS: &[(&str, &[u8])] = &[
    ("Illumina Universal Adapter", b"AGATCGGAAGAG"),
    ("Illumina Small RNA 3' Adapter", b"TGGAATTCTCGG"),
    ("Illumina Small RNA 5' Adapter", b"GATCGTCGGACT"),
    ("Nextera Transposase Sequence", b"CTGTCTCTTATA"),
    ("PolyA", b"AAAAAAAAAAAA"),
    ("PolyG", b"GGGGGGGGGGGG"),
];

pub struct AdapterContent {
    /// Per adapter, per position: count of reads whose first match index
    /// is exactly that position (prefix-summed in `write_data`/status).
    first_hit: Vec<Vec<u64>>,
    total: u64,
    max_len: usize,
}

impl AdapterContent {
    #[must_use]
    pub fn new() -> Self {
        Self {
            first_hit: vec![Vec::new(); ADAPTERS.len()],
            total: 0,
            max_len: 0,
        }
    }

    /// Per-adapter cumulative positive fraction at each position.
    fn cumulative(&self) -> Vec<Vec<f64>> {
        let mut out = Vec::with_capacity(ADAPTERS.len());
        for hits in &self.first_hit {
            let mut row = Vec::with_capacity(self.max_len);
            let mut run: u64 = 0;
            for pos in 0..self.max_len {
                run += hits.get(pos).copied().unwrap_or(0);
                row.push(if self.total == 0 {
                    0.0
                } else {
                    run as f64 / self.total as f64
                });
            }
            out.push(row);
        }
        out
    }
}

impl Default for AdapterContent {
    fn default() -> Self {
        Self::new()
    }
}

impl QcModule for AdapterContent {
    fn name(&self) -> &'static str {
        "Adapter Content"
    }

    fn process(&mut self, rec: &Record) {
        self.total += 1;
        self.max_len = self.max_len.max(rec.seq.len());
        for (ai, &(_, kmer)) in ADAPTERS.iter().enumerate() {
            if rec.seq.len() < kmer.len() {
                continue;
            }
            if let Some(idx) = rec
                .seq
                .windows(kmer.len())
                .position(|w| w.eq_ignore_ascii_case(kmer))
            {
                if idx >= self.first_hit[ai].len() {
                    self.first_hit[ai].resize(idx + 1, 0);
                }
                self.first_hit[ai][idx] += 1;
            }
        }
    }

    fn finalize(&mut self) {}

    fn status(&self) -> ModuleStatus {
        let cum = self.cumulative();
        let mut worst = 0.0_f64;
        for row in &cum {
            for &v in row {
                worst = worst.max(v);
            }
        }
        if worst > 0.10 {
            ModuleStatus::Fail
        } else if worst > 0.05 {
            ModuleStatus::Warn
        } else {
            ModuleStatus::Pass
        }
    }

    fn write_data(&self, out: &mut String) {
        out.push_str("#Position");
        for &(name, _) in ADAPTERS {
            let _ = write!(out, "\t{name}");
        }
        out.push('\n');
        let cum = self.cumulative();
        for pos in 0..self.max_len {
            let _ = write!(out, "{}", pos + 1);
            for row in &cum {
                let _ = write!(out, "\t{:.6}", row[pos] * 100.0);
            }
            out.push('\n');
        }
    }
}
