// u64 counts → f64 only at the final percentage stage.
#![allow(clippy::cast_precision_loss)]

use std::cmp::Reverse;
use std::collections::HashMap;
use std::fmt::Write as _;

use super::{ModuleStatus, QcModule, Record};

/// `FastQC` "Overrepresented sequences" — sequences that are >0.1% of the
/// total are listed; WARN if any sequence exceeds 0.1%, FAIL if any
/// exceeds 1% (clean-room `FastQC` contract). Tracking mirrors the
/// duplication module: only sequences first seen within the first 100 000
/// reads are kept, reads >75 bp keyed by their first 50 bp.
///
/// Each hit is labelled with the best match in `FastQC`'s public
/// `contaminant_list` (a hit must be ≥20 bp with ≤1 mismatch). The full
/// 209-entry public list is embedded before the black-box compat gate; the
/// pass/warn/fail decision is purely count-based and is exact regardless of
/// the label table.
pub struct OverrepresentedSeqs {
    counts: HashMap<Vec<u8>, u64>,
    seen_reads: u64,
    total: u64,
}

const TRACK_LIMIT: u64 = 100_000;
const KEY_TRUNC_OVER: usize = 75;
const KEY_LEN: usize = 50;

/// (name, sequence) — clean-room from `FastQC`'s public `contaminant_list.txt`
/// / `adapter_list.txt` (public config data, not GPL source). Seeded with
/// the common adapters + representative contaminants; completed to the full
/// public list before the compat gate.
const CONTAMINANTS: &[(&str, &[u8])] = &[
    ("Illumina Universal Adapter", b"AGATCGGAAGAG"),
    ("Illumina Small RNA 3' Adapter", b"TGGAATTCTCGG"),
    ("Illumina Small RNA 5' Adapter", b"GATCGTCGGACT"),
    ("Nextera Transposase Sequence", b"CTGTCTCTTATA"),
    ("PolyA", b"AAAAAAAAAAAA"),
    ("PolyG", b"GGGGGGGGGGGG"),
    (
        "TruSeq Universal Adapter",
        b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    ),
    (
        "Illumina Paired End PCR Primer 1",
        b"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    ),
];

fn matches_contaminant(seq: &[u8]) -> Option<&'static str> {
    // FastQC rule: slide the shorter sequence along the longer; a hit is a
    // window with ≤1 mismatch. The overlap must be ≥20 bp, except a
    // contaminant shorter than 20 bp (the 12 bp adapters) must match in
    // full. `short` is whichever of (read, contaminant) is shorter, so its
    // whole length is the overlap.
    for &(name, cont) in CONTAMINANTS {
        let (short, long) = if seq.len() <= cont.len() {
            (seq, cont)
        } else {
            (cont, seq)
        };
        let required = cont.len().min(20);
        if short.is_empty() || short.len() < required {
            continue;
        }
        for start in 0..=(long.len() - short.len()) {
            let window = &long[start..start + short.len()];
            let mismatches = window.iter().zip(short).filter(|(a, b)| a != b).count();
            if mismatches <= 1 {
                return Some(name);
            }
        }
    }
    None
}

impl OverrepresentedSeqs {
    #[must_use]
    pub fn new() -> Self {
        Self {
            counts: HashMap::new(),
            seen_reads: 0,
            total: 0,
        }
    }

    fn key(seq: &[u8]) -> &[u8] {
        if seq.len() > KEY_TRUNC_OVER {
            &seq[..KEY_LEN]
        } else {
            seq
        }
    }

    fn worst_fraction(&self) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        let max = self.counts.values().copied().max().unwrap_or(0);
        max as f64 / self.total as f64
    }
}

impl Default for OverrepresentedSeqs {
    fn default() -> Self {
        Self::new()
    }
}

impl QcModule for OverrepresentedSeqs {
    fn name(&self) -> &'static str {
        "Overrepresented sequences"
    }

    fn process(&mut self, rec: &Record) {
        self.seen_reads += 1;
        self.total += 1;
        let key = Self::key(rec.seq);
        if let Some(c) = self.counts.get_mut(key) {
            *c += 1;
        } else if self.seen_reads <= TRACK_LIMIT {
            self.counts.insert(key.to_vec(), 1);
        }
    }

    fn finalize(&mut self) {}

    fn status(&self) -> ModuleStatus {
        let f = self.worst_fraction();
        if f > 0.01 {
            ModuleStatus::Fail
        } else if f > 0.001 {
            ModuleStatus::Warn
        } else {
            ModuleStatus::Pass
        }
    }

    fn write_data(&self, out: &mut String) {
        out.push_str("#Sequence\tCount\tPercentage\tPossible Source\n");
        if self.total == 0 {
            return;
        }
        let mut rows: Vec<(&Vec<u8>, u64)> = self
            .counts
            .iter()
            .filter(|(_, c)| **c as f64 / self.total as f64 > 0.001)
            .map(|(s, c)| (s, *c))
            .collect();
        rows.sort_by_key(|r| Reverse(r.1));
        for (seq, count) in rows {
            let pct = count as f64 / self.total as f64 * 100.0;
            let source = matches_contaminant(seq).unwrap_or("No Hit");
            let _ = writeln!(
                out,
                "{}\t{count}\t{pct:.6}\t{source}",
                String::from_utf8_lossy(seq)
            );
        }
    }
}
