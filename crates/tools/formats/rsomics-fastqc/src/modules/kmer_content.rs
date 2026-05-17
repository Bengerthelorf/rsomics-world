// u64 counts → f64 only at the final ratio / p-value stage.
#![allow(clippy::cast_precision_loss)]

use std::collections::HashMap;
use std::fmt::Write as _;

use super::{ModuleStatus, QcModule, Record};

/// `FastQC` "Kmer Content" — 7-mers over a 2% read sample (reads truncated
/// to 500 bp). For each k-mer an exact binomial upper-tail test on its most
/// enriched position asks whether it occurs there more than an even spread
/// over positions predicts. WARN if any k-mer's p < 0.01, FAIL if p < 1e-5;
/// the top 6 most biased k-mers are reported (clean-room `FastQC` contract).
///
/// The position support is the widest position span seen in the sample
/// (`n_positions`); for fixed-length reads this is exact. With heavily
/// variable read lengths a short-read k-mer's expectation is taken over the
/// full span, which is conservative — `FastQC`'s exact per-length support
/// model is not in its public documentation, and uniform-length data (the
/// dominant case) is exact.
pub struct KmerContent {
    /// k-mer (2-bit packed, 7 bases) → per-position observed counts.
    counts: HashMap<u32, Vec<u64>>,
    read_idx: u64,
    n_positions: usize,
    stats_cache: Vec<KmerStat>,
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

/// Log-gamma via the Lanczos approximation (g=5, 6 terms) — standard
/// public-domain coefficients; accurate to ~1e-10 for the a,b used here.
fn gammaln(x: f64) -> f64 {
    const C: [f64; 6] = [
        76.180_091_729_471_46,
        -86.505_320_329_416_77,
        24.014_098_240_830_91,
        -1.231_739_572_450_155,
        0.001_208_650_973_866_179,
        -0.000_005_395_239_384_953,
    ];
    let mut ser = 1.000_000_000_190_015;
    let mut y = x;
    for c in C {
        y += 1.0;
        ser += c / y;
    }
    let tmp = x + 5.5 - (x + 0.5) * (x + 5.5).ln();
    -tmp + (2.506_628_274_631_000_5 * ser / x).ln()
}

/// Continued fraction for the incomplete beta (modified Lentz). The
/// recurrence is the canonical textbook method, written from the math;
/// the single-letter names are its standard mathematical notation
/// (a, b, x, and the Lentz state c, d, h) — verbose names would obscure
/// the algorithm.
#[allow(clippy::many_single_char_names)]
fn betacf(a: f64, b: f64, x: f64) -> f64 {
    let tiny = 1e-30;
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < tiny {
        d = tiny;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..200 {
        let m = f64::from(m);
        let m2 = 2.0 * m;
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < tiny {
            d = tiny;
        }
        c = 1.0 + aa / c;
        if c.abs() < tiny {
            c = tiny;
        }
        d = 1.0 / d;
        h *= d * c;
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < tiny {
            d = tiny;
        }
        c = 1.0 + aa / c;
        if c.abs() < tiny {
            c = tiny;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-12 {
            break;
        }
    }
    h
}

/// Regularised incomplete beta `I_x(a, b)`.
fn betai(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let bt = (gammaln(a + b) - gammaln(a) - gammaln(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        bt * betacf(a, b, x) / a
    } else {
        1.0 - bt * betacf(b, a, 1.0 - x) / b
    }
}

/// Exact binomial upper tail `P(X >= k)`, `X ~ Binomial(n, p)`, via the
/// identity `P(X >= k) = I_p(k, n-k+1)` (k >= 1). `FastQC` uses the exact
/// binomial — a normal approximation is invalid where the per-position
/// expectation is small (`n*p << 5`), the regime that decides WARN/FAIL.
fn binom_sf(k: u64, n: u64, p: f64) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if k > n {
        return 0.0;
    }
    betai(k as f64, (n - k + 1) as f64, p)
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
            stats_cache: Vec::new(),
        }
    }

    fn compute_stats(&self) -> Vec<KmerStat> {
        let mut v = Vec::with_capacity(self.counts.len());
        for (&code, per_pos) in &self.counts {
            let total: u64 = per_pos.iter().sum();
            if total == 0 || self.n_positions == 0 {
                continue;
            }
            let p0 = 1.0 / self.n_positions as f64;
            let mean = total as f64 * p0;
            let mut max_pos = 0usize;
            let mut max_obs = 0u64;
            for (i, &c) in per_pos.iter().enumerate() {
                if c > max_obs {
                    max_obs = c;
                    max_pos = i;
                }
            }
            v.push(KmerStat {
                seq: unpack(code),
                count: total,
                p: binom_sf(max_obs, total, p0),
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

    fn finalize(&mut self) {
        self.stats_cache = self.compute_stats();
    }

    fn status(&self) -> ModuleStatus {
        let mut worst_p = 1.0_f64;
        for s in &self.stats_cache {
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
        for s in self.stats_cache.iter().take(6) {
            let _ = writeln!(
                out,
                "{}\t{}\t{:.6E}\t{:.6}\t{}",
                s.seq, s.count, s.p, s.obs_exp_max, s.max_pos
            );
        }
    }
}
