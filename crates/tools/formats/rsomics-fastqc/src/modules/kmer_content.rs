// u64 counts → f64 only at the final ratio / p-value stage.
#![allow(clippy::cast_precision_loss)]

use std::collections::HashMap;
use std::fmt::Write as _;

use super::{ModuleStatus, QcModule, Record};

/// `FastQC` "Kmer Content" — 7-mers over a 2% read sample (reads truncated to
/// 500 bp). For each k-mer a binomial test looks for a position with more
/// occurrences than an even spread predicts. WARN if any k-mer's most
/// enriched position has binomial p < 0.01, FAIL if p < 1e-5; the top 6
/// most biased k-mers are reported (clean-room `FastQC` contract).
///
/// `FastQC` uses an exact binomial; this uses the normal approximation with
/// continuity correction (n is large) — exact p-values are calibrated
/// against the `FastQC` binary at the black-box compat step, the pass/warn/
/// fail thresholds are implemented per the documented contract.
pub struct KmerContent {
    /// k-mer (2-bit packed, 7 bases) → per-position observed counts.
    counts: HashMap<u32, Vec<u64>>,
    read_idx: u64,
    n_positions: usize,
}

const K: usize = 7;
const SAMPLE_EVERY: u64 = 50; // 2% of reads
const MAX_LEN: usize = 500;

fn base_2bit(b: u8) -> Option<u32> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn unpack(mut code: u32) -> String {
    let mut s = [0u8; K];
    for i in (0..K).rev() {
        s[i] = b"ACGT"[(code & 3) as usize];
        code >>= 2;
    }
    String::from_utf8_lossy(&s).into_owned()
}

/// Upper-tail standard-normal probability via `erfc` (Abramowitz & Stegun
/// 7.1.26) — the binomial normal approximation's survival function.
fn norm_sf(z: f64) -> f64 {
    let x = z / std::f64::consts::SQRT_2;
    let t = 1.0 / (1.0 + 0.327_591_1 * x.abs());
    let y = 1.0
        - (((((1.061_405_429 * t - 1.453_152_027) * t) + 1.421_413_741) * t - 0.284_496_736) * t
            + 0.254_829_592)
            * t
            * (-x * x).exp();
    let erf = if x >= 0.0 { y } else { -y };
    0.5 * (1.0 - erf)
}

struct KmerStat {
    seq: String,
    count: u64,
    p: f64,
    obs_exp_max: f64,
    max_pos: usize,
}

impl KmerContent {
    #[must_use]
    pub fn new() -> Self {
        Self {
            counts: HashMap::new(),
            read_idx: 0,
            n_positions: 0,
        }
    }

    fn stats(&self) -> Vec<KmerStat> {
        let mut v = Vec::with_capacity(self.counts.len());
        for (&code, per_pos) in &self.counts {
            let total: u64 = per_pos.iter().sum();
            if total == 0 || self.n_positions == 0 {
                continue;
            }
            let p0 = 1.0 / self.n_positions as f64;
            let mean = total as f64 * p0;
            let sd = (total as f64 * p0 * (1.0 - p0)).sqrt().max(1e-9);
            let mut max_pos = 0usize;
            let mut max_obs = 0u64;
            for (i, &c) in per_pos.iter().enumerate() {
                if c > max_obs {
                    max_obs = c;
                    max_pos = i;
                }
            }
            let z = (max_obs as f64 - mean - 0.5) / sd;
            v.push(KmerStat {
                seq: unpack(code),
                count: total,
                p: norm_sf(z),
                obs_exp_max: if mean > 0.0 {
                    max_obs as f64 / mean
                } else {
                    0.0
                },
                max_pos: max_pos + 1,
            });
        }
        v.sort_by(|a, b| b.obs_exp_max.total_cmp(&a.obs_exp_max));
        v
    }
}

impl Default for KmerContent {
    fn default() -> Self {
        Self::new()
    }
}

impl QcModule for KmerContent {
    fn name(&self) -> &'static str {
        "Kmer Content"
    }

    fn process(&mut self, rec: &Record) {
        self.read_idx += 1;
        if self.read_idx % SAMPLE_EVERY != 1 {
            return;
        }
        let seq = &rec.seq[..rec.seq.len().min(MAX_LEN)];
        if seq.len() < K {
            return;
        }
        let last_start = seq.len() - K;
        self.n_positions = self.n_positions.max(last_start + 1);
        let npos = self.n_positions;
        for start in 0..=last_start {
            let mut code = 0u32;
            let mut ok = true;
            for &b in &seq[start..start + K] {
                if let Some(c) = base_2bit(b) {
                    code = (code << 2) | c;
                } else {
                    ok = false;
                    break;
                }
            }
            if !ok {
                continue;
            }
            let row = self.counts.entry(code).or_insert_with(|| vec![0; npos]);
            if start >= row.len() {
                row.resize(start + 1, 0);
            }
            row[start] += 1;
        }
    }

    fn finalize(&mut self) {}

    fn status(&self) -> ModuleStatus {
        let mut worst_p = 1.0_f64;
        for s in self.stats() {
            worst_p = worst_p.min(s.p);
        }
        if worst_p < 1e-5 {
            ModuleStatus::Fail
        } else if worst_p < 0.01 {
            ModuleStatus::Warn
        } else {
            ModuleStatus::Pass
        }
    }

    fn write_data(&self, out: &mut String) {
        out.push_str("#Sequence\tCount\tPValue\tObs/Exp Max\tMax Obs/Exp Position\n");
        for s in self.stats().into_iter().take(6) {
            let _ = writeln!(
                out,
                "{}\t{}\t{:.6E}\t{:.6}\t{}",
                s.seq, s.count, s.p, s.obs_exp_max, s.max_pos
            );
        }
    }
}
