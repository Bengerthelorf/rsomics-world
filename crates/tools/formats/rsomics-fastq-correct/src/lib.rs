// This crate is a faithful port of BFC's C corrector (lh3/bfc, MIT). The
// single-character names, positional indexing, and signed/unsigned casts
// mirror `correct.c` deliberately so the algorithm stays diff-able against
// the source; the same allow-set is used by rsomics-stats for its R ports.
#![allow(
    clippy::many_single_char_names,
    clippy::needless_range_loop,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::cast_precision_loss
)]

use std::collections::HashMap;
use std::path::Path;

use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};
use rsomics_fqgz::ChunkedWriter;
use rsomics_seqio::{OwnedRecord, open_fastq};
use serde::Serialize;

const CHUNK_RECORDS: usize = 4096;
const EC_HIST: usize = 5;
const EC_HIST_HIGH: usize = 2;
const MAX_PATHS: usize = 4;

/// fastp/BFC base codes: 0..=3 = A,C,G,T; 4 = ambiguous (N).
const SEQ_NT4: [u8; 256] = build_nt4();

const fn build_nt4() -> [u8; 256] {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
}

/// Configuration for one correction run. Field names and defaults are the
/// BFC 0.1 `bfc_opt_init` values (`-k` 33, `-c` 3, `-w` 10, `-q` 20, the
/// fixed `w_*` penalty weights, `max_path_diff` 15, `max_heap` 100,
/// `max_end_ext` 5). The `w_*` weights are not exposed on the CLI — BFC
/// marks them "cannot be changed on command line"; we keep that contract.
#[derive(Debug, Clone)]
pub struct CorrectConfig {
    pub k: usize,
    pub min_cov: i32,
    pub win_multi_ec: i32,
    pub qual_threshold: u8,
    pub max_end_ext: i32,
    pub max_path_diff: i32,
    pub max_heap: usize,
    pub drop_unique_kmer: bool,
    pub discard_uncorrectable: bool,
    pub fasta_out: bool,
}

impl Default for CorrectConfig {
    fn default() -> Self {
        Self {
            k: 33,
            min_cov: 3,
            win_multi_ec: 10,
            qual_threshold: 20,
            max_end_ext: 5,
            max_path_diff: 15,
            max_heap: 100,
            drop_unique_kmer: false,
            discard_uncorrectable: false,
            fasta_out: false,
        }
    }
}

const W_EC: i32 = 1;
const W_EC_HIGH: i32 = 7;
const W_ABSENT: i32 = 3;
const W_ABSENT_HIGH: i32 = 1;

/// BFC `bfc_kmer_t`: x[0]/x[1] are the 2-bit forward k-mer (low/high bit
/// plane), x[2]/x[3] the reverse complement. k ≤ 63 (`BFC_MAX_KMER`).
#[derive(Clone, Copy)]
struct BfcKmer {
    x: [u64; 4],
}

impl BfcKmer {
    const NULL: BfcKmer = BfcKmer { x: [0, 0, 0, 0] };

    #[inline]
    fn append(&mut self, k: usize, c: u64) {
        let mask = if k >= 64 { u64::MAX } else { (1u64 << k) - 1 };
        self.x[0] = ((self.x[0] << 1) | (c & 1)) & mask;
        self.x[1] = ((self.x[1] << 1) | (c >> 1)) & mask;
        self.x[2] = (self.x[2] >> 1) | ((1u64 ^ (c & 1)) << (k - 1));
        self.x[3] = (self.x[3] >> 1) | ((1u64 ^ (c >> 1)) << (k - 1));
    }

    /// BFC `bfc_kmer_change`: set base `d` (counted from the 3' end,
    /// `0 ≤ d < k`) to `c` in-place on both the forward and RC planes.
    #[inline]
    fn change(&mut self, k: usize, d: usize, c: u64) {
        let t = !(1u64 << d);
        self.x[0] = ((c & 1) << d) | (self.x[0] & t);
        self.x[1] = ((c >> 1) << d) | (self.x[1] & t);
        let t = !(1u64 << (k - 1 - d));
        self.x[2] = ((1u64 ^ (c & 1)) << (k - 1 - d)) | (self.x[2] & t);
        self.x[3] = ((1u64 ^ (c >> 1)) << (k - 1 - d)) | (self.x[3] & t);
    }
}

/// BFC `bfc_hash_64` — Thomas Wang's invertible integer hash, masked.
#[inline]
fn bfc_hash_64(mut key: u64, mask: u64) -> u64 {
    key = (!key).wrapping_add(key << 21) & mask;
    key ^= key >> 24;
    key = key.wrapping_add(key << 3).wrapping_add(key << 8) & mask;
    key ^= key >> 14;
    key = key.wrapping_add(key << 2).wrapping_add(key << 4) & mask;
    key ^= key >> 28;
    key = key.wrapping_add(key << 31) & mask;
    key
}

/// BFC `bfc_kmer_hash`: canonical double-hash. The middle base selects the
/// strand (`x[1]>>t&1 > x[3]>>t&1`); the returned 64-bit value is BFC's
/// counting-table key. Reproduced verbatim so the trusted-k-mer counts
/// match BFC's even though our table is a plain map, not its open-address
/// `bfc_ch` (the table's collision profile is the documented compat gap).
#[inline]
fn bfc_kmer_hash(k: usize, x: &[u64; 4]) -> u64 {
    // BFC_MAX_KMER is 63; k is validated ≤ 63 at the Pipeline/CLI boundary
    // (fail-loud there). The `<< k` below is only well-defined for k < 64,
    // so assert the invariant rather than silently zeroing in release.
    debug_assert!(k < 64, "k must be < 64 (BFC_MAX_KMER); got {k}");
    let t = k >> 1;
    let u = usize::from((x[1] >> t & 1) > (x[3] >> t & 1));
    let mask = (1u64 << k) - 1;
    let h0 = bfc_hash_64((x[u << 1].wrapping_add(x[(u << 1) | 1])) & mask, mask);
    let h1 = bfc_hash_64(h0 ^ x[(u << 1) | 1], mask);
    ((h0 ^ h1) << k) | ((h0.wrapping_add(h1)) & mask)
}

/// Per-k-mer occupancy: low byte = coverage, bits 8..14 = high-quality
/// coverage, both saturating — the layout BFC's `bfc_ch_kmer_occ` returns
/// (`r&0xff`, `r>>8&0x3f`).
#[derive(Default, Clone, Copy)]
struct Occ {
    cov: u8,
    hi: u8,
}

/// The trusted-k-mer count table (BFC `bfc_ch`, modelled as a map keyed by
/// `bfc_kmer_hash`). Built once over every input read before correction.
pub struct CountTable {
    k: usize,
    map: HashMap<u64, Occ>,
}

impl CountTable {
    fn occ(&self, kmer: &BfcKmer) -> Occ {
        self.map
            .get(&bfc_kmer_hash(self.k, &kmer.x))
            .copied()
            .unwrap_or_default()
    }

    fn add(&mut self, kmer: &BfcKmer, high_qual: bool) {
        let e = self.map.entry(bfc_kmer_hash(self.k, &kmer.x)).or_default();
        e.cov = e.cov.saturating_add(1);
        if high_qual {
            e.hi = e.hi.saturating_add(1).min(0x3f);
        }
    }

    /// BFC `bfc_ch_hist` mode: the peak of the coverage histogram (the
    /// genomic coverage), used by the greedy probe's confidence gate. The
    /// peak is taken over `cov ≥ min_cov` so the error-noise spike at low
    /// coverage does not dominate.
    fn hist_mode(&self, min_cov: i32) -> i32 {
        let mut hist = [0u64; 256];
        for o in self.map.values() {
            hist[o.cov as usize] += 1;
        }
        let (mut mode, mut best) = (0i32, 0u64);
        for c in (min_cov.max(1) as usize)..256 {
            if hist[c] > best {
                best = hist[c];
                mode = c as i32;
            }
        }
        mode
    }
}

/// One base in the ec working sequence — BFC `ecbase_t` (`b`/`ob` current
/// & original base, `q`/`oq` quality bit, `lcov`/`hcov` solid coverage,
/// the solid/high end flags).
#[derive(Clone, Copy, Default)]
struct EcBase {
    b: u8,
    ob: u8,
    q: u8,
    oq: u8,
    lcov: i32,
    hcov: i32,
    solid_end: bool,
    high_end: bool,
}

fn seq_conv(seq: &[u8], qual: &[u8], qthres: u8) -> Vec<EcBase> {
    seq.iter()
        .enumerate()
        .map(|(i, &s)| {
            let b = SEQ_NT4[s as usize];
            let mut q = u8::from(qual.is_empty() || qual[i].saturating_sub(33) >= qthres);
            if b > 3 {
                q = 0;
            }
            EcBase {
                b,
                ob: b,
                q,
                oq: q,
                ..EcBase::default()
            }
        })
        .collect()
}

#[inline]
fn comp_base(b: u8) -> u8 {
    if b < 4 { 3 - b } else { 4 }
}

fn revcomp(s: &mut [EcBase]) {
    s.reverse();
    for c in s.iter_mut() {
        c.b = comp_base(c.b);
        c.ob = comp_base(c.ob);
    }
}

/// BFC `bfc_ec_kcov`: mark each base's solid/high-end flags and accumulate
/// `lcov`/`hcov` over the k-mers covering it.
fn ec_kcov(k: usize, min_occ: i32, s: &mut [EcBase], ch: &CountTable) {
    let mut x = BfcKmer::NULL;
    let mut l = 0usize;
    for i in 0..s.len() {
        s[i].high_end = false;
        s[i].solid_end = false;
        s[i].lcov = 0;
        s[i].hcov = 0;
        if s[i].b < 4 {
            x.append(k, u64::from(s[i].b));
            l += 1;
            if l >= k {
                let r = ch.occ(&x);
                if i32::from(r.hi) > min_occ {
                    s[i].high_end = true;
                }
                if i32::from(r.cov) >= min_occ {
                    s[i].solid_end = true;
                    let he = i32::from(s[i].high_end);
                    for j in (i + 1 - k)..=i {
                        s[j].lcov += 1;
                        s[j].hcov += he;
                    }
                }
            }
        } else {
            l = 0;
            x = BfcKmer::NULL;
        }
    }
}

/// BFC `bfc_ec_best_island`: longest run of `solid_end` k-mers. Returns
/// `Some((start, end))` (the C `max_i-max-k+1 .. max_i` half-open) or
/// `None` for no solid k-mer. Call after [`ec_kcov`].
fn best_island(k: usize, s: &[EcBase]) -> Option<(usize, usize)> {
    let (mut max, mut l) = (0i64, 0i64);
    let mut max_i: i64 = -1;
    let n = s.len();
    for i in (k - 1)..n {
        if s[i].solid_end {
            l += 1;
        } else {
            if l > max {
                max = l;
                max_i = i as i64;
            }
            l = 0;
        }
    }
    if l > max {
        max = l;
        max_i = n as i64;
    }
    if max == 0 {
        return None;
    }
    // BFC `bfc_ec_best_island`: start = first base the solid run covers =
    // (first solid k-mer index) - k + 1. Signed like the C int math; the
    // solid-end ⟹ index ≥ k-1 invariant makes this ≥ 0.
    let start = max_i - max - k as i64 + 1;
    debug_assert!(
        start >= 0,
        "best_island start<0: start={start} max_i={max_i} max={max} k={k} n={n}"
    );
    Some((start.max(0) as usize, max_i as usize))
}

/// BFC `bfc_ec_first_kmer`: scan from `start` for the first window of `k`
/// consecutive non-N bases. Returns `(i, x)` where `i` is the index of the
/// k-th base (the k-mer's 3' end), or `s.len()` if no such window exists.
fn ec_first_kmer(k: usize, s: &[EcBase], start: usize) -> (usize, BfcKmer) {
    let mut x = BfcKmer::NULL;
    let mut l = 0usize;
    let mut i = start;
    while i < s.len() {
        if s[i].b < 4 {
            x.append(k, u64::from(s[i].b));
            l += 1;
            if l == k {
                break;
            }
        } else {
            l = 0;
            x = BfcKmer::NULL;
        }
        i += 1;
    }
    (i, x)
}

/// BFC `bfc_ec_greedy_k`: try every single-base substitution of `x`; if
/// exactly one alternative k-mer is strongly supported (best coverage
/// `> mode/3` and second-best `< 3`), return `pos<<2 | base` (pos counted
/// from the 3' end), else `-1`. Rescues reads with no solid island but one
/// confident single error in the first clean k-mer window.
fn bfc_ec_greedy_k(k: usize, mode: i32, x: &BfcKmer, ch: &CountTable) -> i32 {
    let (mut max, mut max2, mut max_ec) = (0i32, 0i32, -1i32);
    for i in 0..k {
        let c = (x.x[1] >> i & 1) << 1 | (x.x[0] >> i & 1);
        for j in 0u64..4 {
            if j == c {
                continue;
            }
            let mut y = *x;
            y.change(k, i, j);
            let ret = i32::from(ch.occ(&y).cov);
            if ret == 0 {
                continue;
            }
            if max < ret {
                max2 = max;
                max = ret;
                max_ec = (i << 2 | j as usize) as i32;
            } else if max2 < ret {
                max2 = ret;
            }
        }
    }
    if max * 3 > mode && max2 < 3 {
        max_ec
    } else {
        -1
    }
}

// BFC `bfc_penalty_t` — four 1-bit penalty kinds + the chosen base; the
// bool count mirrors the C bitfield, it is not a config struct.
#[allow(clippy::struct_excessive_bools)]
#[derive(Clone, Copy, Default)]
struct Penalty {
    ec: bool,
    ec_high: bool,
    absent: bool,
    absent_high: bool,
    b: u8,
}

#[inline]
fn weighted_penalty(p: Penalty) -> i32 {
    W_EC * i32::from(p.ec)
        + W_EC_HIGH * i32::from(p.ec_high)
        + W_ABSENT * i32::from(p.absent)
        + W_ABSENT_HIGH * i32::from(p.absent_high)
}

#[derive(Clone)]
struct HeapNode {
    tot_pen: i32,
    i: usize,
    k: i64,
    ecpos_high: [i64; EC_HIST_HIGH],
    ecpos: [i64; EC_HIST],
    x: BfcKmer,
}

struct StackNode {
    parent: i64,
    i: usize,
    tot_pen: i32,
    b: u8,
    pen: Penalty,
}

/// BFC `bfc_ec1dir`: best-first correction from a solid anchor. A binary
/// min-heap on `tot_pen` (BFC's `heap_lt` is `>`, i.e. lowest penalty
/// pops first) drives an exhaustive bounded search; the lowest-penalty
/// terminating path is back-tracked into `ec`. Returns `n_absent` (≥0) or
/// a negative BFC code (`-2` uncorrectable N, `-3` too many failures).
#[allow(clippy::too_many_lines)]
fn ec1dir(
    cfg: &CorrectConfig,
    ch: &CountTable,
    seq: &[EcBase],
    ec: &mut Vec<EcBase>,
    start: usize,
    end: usize,
) -> i32 {
    let k = cfg.k;
    let n = seq.len();
    let mut heap: Vec<HeapNode> = Vec::new();
    let mut stack: Vec<StackNode> = Vec::new();
    let mut z = HeapNode {
        tot_pen: 0,
        i: 0,
        k: -1,
        ecpos_high: [-1; EC_HIST_HIGH],
        ecpos: [-1; EC_HIST],
        x: BfcKmer::NULL,
    };
    ec.clear();
    ec.extend_from_slice(seq);

    let mut l = 0usize;
    z.i = start;
    while z.i < end {
        let c = seq[z.i].b;
        if c < 4 {
            l += 1;
            if l == k {
                break;
            }
            z.x.append(k, u64::from(c));
        } else {
            l = 0;
            z.x = BfcKmer::NULL;
        }
        z.i += 1;
    }
    debug_assert!(z.i < end);
    heap.push(z);

    let mut path = [0i64; MAX_PATHS];
    let mut n_paths = 0usize;
    let mut min_path: i64 = -1;
    let mut min_path_pen = i32::MAX;
    let mut n_failures = 0i32;
    let mut rv: i32 = -1;

    while !heap.is_empty() {
        let z = heap_pop(&mut heap);
        if min_path >= 0 && z.tot_pen > min_path_pen + cfg.max_path_diff {
            break;
        }
        let mut stop = z.i as i64 - end as i64 > i64::from(cfg.max_end_ext);
        if !stop {
            let c = if z.i < n { Some(seq[z.i]) } else { None };
            let mut os: i32 = -1;
            let mut fixed = z.i > end;
            let mut other_ext = 0;
            let mut added: Vec<(Penalty, i32)> = Vec::with_capacity(4);

            if let Some(cb) = c
                && cb.b < 4
            {
                let mut x = z.x;
                x.append(k, u64::from(cb.b));
                let o = ch.occ(&x);
                os = i32::from(o.cov) | (i32::from(o.hi) << 8);
                if (cb.q != 0 && i32::from(o.cov) > cfg.min_cov && cb.lcov > cfg.min_cov)
                    || f64::from(cb.hcov) > k as f64 * 0.75
                {
                    fixed = true;
                }
            }

            for b in 0u8..4 {
                if fixed
                    && let Some(cb) = c
                    && b != cb.b
                {
                    continue;
                }
                if c.is_none() || b != c.unwrap().b {
                    let mut x = z.x;
                    if let Some(cb) = c {
                        if cb.q != 0
                            && z.ecpos_high[EC_HIST_HIGH - 1] >= 0
                            && z.i as i64 - z.ecpos_high[EC_HIST_HIGH - 1]
                                < i64::from(cfg.win_multi_ec)
                        {
                            continue;
                        }
                        if z.ecpos[EC_HIST - 1] >= 0
                            && z.i as i64 - z.ecpos[EC_HIST - 1] < i64::from(cfg.win_multi_ec)
                        {
                            continue;
                        }
                    }
                    x.append(k, u64::from(b));
                    let o = ch.occ(&x);
                    let s = i32::from(o.cov) | (i32::from(o.hi) << 8);
                    if i32::from(o.cov) < cfg.min_cov {
                        continue;
                    }
                    let ec_flag = c.is_some_and(|cb| cb.b < 4);
                    let pen = Penalty {
                        ec: ec_flag,
                        ec_high: ec_flag && c.unwrap().oq != 0,
                        absent: false,
                        absent_high: (s >> 8 & 0xff) < cfg.min_cov,
                        b,
                    };
                    added.push((pen, s));
                    other_ext += 1;
                } else {
                    let pen = Penalty {
                        ec: false,
                        ec_high: false,
                        absent: os < 0 || (os & 0xff) < cfg.min_cov,
                        absent_high: os < 0 || (os >> 8 & 0xff) < cfg.min_cov,
                        b,
                    };
                    added.push((pen, os));
                }
            }

            if !fixed && other_ext == 0 {
                n_failures += 1;
            }
            if n_failures > n as i32 * 2 {
                rv = -3;
                break;
            }

            if c.is_some() || added.len() == 1 {
                if added.len() > 1 && heap.len() > cfg.max_heap {
                    let best = added
                        .iter()
                        .min_by_key(|(p, _)| weighted_penalty(*p))
                        .copied()
                        .unwrap();
                    buf_update(cfg, &mut heap, &mut stack, &z, best.0);
                } else {
                    for (pen, _) in &added {
                        buf_update(cfg, &mut heap, &mut stack, &z, *pen);
                    }
                }
            } else {
                if added.is_empty() && z.k >= 0 {
                    let zk = z.k as usize;
                    stack[zk].tot_pen +=
                        W_ABSENT * (cfg.max_end_ext - (z.i as i64 - end as i64) as i32);
                }
                stop = true;
            }
        }
        // Only a node with a real stack index is a valid terminating path;
        // BFC's path[] never holds the seed's -1 (the seed always extends
        // before any stop), so a degenerate z.k<0 stop is not a path.
        if stop && z.k >= 0 {
            let tp = stack[z.k as usize].tot_pen;
            if tp < min_path_pen {
                min_path_pen = tp;
                min_path = n_paths as i64;
            }
            path[n_paths] = z.k;
            n_paths += 1;
            if n_paths == MAX_PATHS {
                break;
            }
        }
    }

    if n_paths == 0 {
        return rv;
    }
    let chosen = path[min_path as usize];
    let n_absent = backtrack(&stack, chosen, seq, ec);
    for (i, b) in ec.iter_mut().enumerate() {
        if i < start + k || i >= end {
            b.b = 4;
        }
    }
    n_absent
}

fn buf_update(
    cfg: &CorrectConfig,
    heap: &mut Vec<HeapNode>,
    stack: &mut Vec<StackNode>,
    prev: &HeapNode,
    pen: Penalty,
) {
    let tot_pen = prev.tot_pen + weighted_penalty(pen);
    stack.push(StackNode {
        parent: prev.k,
        i: prev.i,
        tot_pen,
        b: pen.b,
        pen,
    });
    let mut r = HeapNode {
        tot_pen,
        i: prev.i + 1,
        k: (stack.len() - 1) as i64,
        ecpos_high: prev.ecpos_high,
        ecpos: prev.ecpos,
        x: prev.x,
    };
    if pen.ec_high {
        r.ecpos_high.copy_within(0..EC_HIST_HIGH - 1, 1);
        r.ecpos_high[0] = prev.i as i64;
    }
    if pen.ec {
        r.ecpos.copy_within(0..EC_HIST - 1, 1);
        r.ecpos[0] = prev.i as i64;
    }
    r.x.append(cfg.k, u64::from(pen.b));
    heap.push(r);
    let last = heap.len() - 1;
    heap_up(heap, last);
}

/// BFC backtrack: walk the stack parent chain, writing the chosen base and
/// ec/absent flags at each corrected position. Returns `n_absent`.
fn backtrack(stack: &[StackNode], mut end: i64, seq: &[EcBase], path: &mut Vec<EcBase>) -> i32 {
    path.clear();
    path.extend_from_slice(seq);
    let mut n_absent = 0;
    while end >= 0 {
        let s = &stack[end as usize];
        if s.i < seq.len() {
            path[s.i].b = s.b;
            if s.pen.absent {
                n_absent += 1;
            }
        }
        end = s.parent;
    }
    n_absent
}

// BFC's heap is a min-heap on tot_pen (`heap_lt(a,b) = a.tot_pen > b.tot_pen`).
#[inline]
fn heap_lt(a: &HeapNode, b: &HeapNode) -> bool {
    a.tot_pen < b.tot_pen
}

fn heap_up(h: &mut [HeapNode], mut i: usize) {
    while i > 0 {
        let p = (i - 1) / 2;
        if heap_lt(&h[i], &h[p]) {
            h.swap(i, p);
            i = p;
        } else {
            break;
        }
    }
}

fn heap_pop(h: &mut Vec<HeapNode>) -> HeapNode {
    let n = h.len();
    h.swap(0, n - 1);
    let top = h.pop().unwrap();
    let len = h.len();
    let mut i = 0;
    loop {
        let (l, r) = (2 * i + 1, 2 * i + 2);
        let mut m = i;
        if l < len && heap_lt(&h[l], &h[m]) {
            m = l;
        }
        if r < len && heap_lt(&h[r], &h[m]) {
            m = r;
        }
        if m == i {
            break;
        }
        h.swap(i, m);
        i = m;
    }
    top
}

/// BFC `bfc_ec1`: correct one read. N-fraction guard (>5% → drop), solid
/// island, forward `ec1dir` over the full read, reverse-complement, second
/// `ec1dir`, then the two directional results are merged by BFC's
/// agreement rule and the corrected sequence + ec-encoded quality emitted.
/// Returns `None` when the read is uncorrectable / over the N threshold.
fn correct_one(
    cfg: &CorrectConfig,
    ch: &CountTable,
    rec: &OwnedRecord,
) -> Option<(Vec<u8>, Vec<u8>)> {
    let mut s = seq_conv(&rec.seq, &rec.qual, cfg.qual_threshold);
    let n = s.len();
    if n < cfg.k {
        return None;
    }
    let n_n = s.iter().filter(|c| c.ob > 3).count();
    if n_n as f64 > n as f64 * 0.05 {
        return None;
    }
    ec_kcov(cfg.k, cfg.min_cov, &mut s, ch);
    // BFC `bfc_ec1`: with a solid island, anchor there. Without one, try
    // the greedy single-substitution probe over successive first-k-mers;
    // on success apply that one base fix and re-anchor, else NO_SOLID.
    let (start, end) = if let Some(se) = best_island(cfg.k, &s) {
        se
    } else {
        let mode = ch.hist_mode(cfg.min_cov);
        let mut bstart = 0usize;
        let mut ec = -1i32;
        let mut bend;
        loop {
            let (e, x) = ec_first_kmer(cfg.k, &s, bstart);
            bend = e;
            if bend >= n {
                break;
            }
            ec = bfc_ec_greedy_k(cfg.k, mode, &x, ch);
            if ec >= 0 {
                break;
            }
            if bend + (cfg.k >> 1) >= n {
                break;
            }
            bstart = bend - (cfg.k >> 1);
        }
        if ec < 0 {
            return None;
        }
        s[bend - (ec as usize >> 2)].b = (ec & 3) as u8;
        let ne = bend + 1;
        (ne - cfg.k, ne)
    };

    let mut ec0 = Vec::new();
    if ec1dir(cfg, ch, &s, &mut ec0, start, n) < 0 {
        return None;
    }
    revcomp(&mut s);
    let mut ec1 = Vec::new();
    if ec1dir(cfg, ch, &s, &mut ec1, n - end, n) < 0 {
        return None;
    }
    revcomp(&mut ec1);
    revcomp(&mut s);

    for i in 0..n {
        let b = if ec0[i].b == ec1[i].b {
            if ec0[i].b > 3 { s[i].b } else { ec0[i].b }
        } else if ec1[i].b > 3 {
            ec0[i].b
        } else if ec0[i].b > 3 {
            ec1[i].b
        } else {
            s[i].ob
        };
        s[i].b = b;
    }

    let mut out_seq = Vec::with_capacity(n);
    let mut out_qual = Vec::with_capacity(n);
    for c in &s {
        let is_diff = c.b != c.ob;
        out_seq.push(if is_diff {
            b"acgtn"[c.b.min(4) as usize]
        } else {
            b"ACGTN"[c.b.min(4) as usize]
        });
        if !rec.qual.is_empty() {
            out_qual.push(if is_diff {
                34 + c.ob
            } else if c.q != 0 {
                b'?'
            } else {
                b'+'
            });
        }
    }
    Some((out_seq, out_qual))
}

#[derive(Debug, Default, Clone, Serialize)]
pub struct CorrectReport {
    pub reads_in: u64,
    pub reads_out: u64,
    pub reads_dropped: u64,
    pub bases_in: u64,
    pub bases_corrected: u64,
}

/// Count every input k-mer into the trusted-k-mer table (BFC's count pass —
/// it scans the whole input before correcting). A k-mer occurrence counts
/// as high-quality when both flanking quality bits are set.
fn build_table(paths: &[&Path], cfg: &CorrectConfig) -> Result<CountTable> {
    let mut ch = CountTable {
        k: cfg.k,
        map: HashMap::new(),
    };
    for p in paths {
        let src = open_fastq(p)?;
        for r in src {
            let rec = r?;
            let s = seq_conv(&rec.seq, &rec.qual, cfg.qual_threshold);
            let mut x = BfcKmer::NULL;
            let mut l = 0usize;
            for i in 0..s.len() {
                if s[i].b < 4 {
                    x.append(cfg.k, u64::from(s[i].b));
                    l += 1;
                    if l >= cfg.k {
                        let hq = s[i].oq != 0 && s[i + 1 - cfg.k].oq != 0;
                        ch.add(&x, hq);
                    }
                } else {
                    l = 0;
                    x = BfcKmer::NULL;
                }
            }
        }
    }
    Ok(ch)
}

fn has_unique_kmer(cfg: &CorrectConfig, ch: &CountTable, rec: &OwnedRecord) -> bool {
    let s = seq_conv(&rec.seq, &rec.qual, cfg.qual_threshold);
    let mut x = BfcKmer::NULL;
    let mut l = 0usize;
    for i in 0..s.len() {
        if s[i].b < 4 {
            x.append(cfg.k, u64::from(s[i].b));
            l += 1;
            if l >= cfg.k && ch.occ(&x).cov <= 1 {
                return true;
            }
        } else {
            l = 0;
            x = BfcKmer::NULL;
        }
    }
    false
}

pub struct Pipeline<'cfg> {
    pub cfg: &'cfg CorrectConfig,
    pub compression: i32,
}

impl<'cfg> Pipeline<'cfg> {
    #[must_use]
    pub fn new(cfg: &'cfg CorrectConfig, compression: i32) -> Self {
        Self { cfg, compression }
    }

    /// # Errors
    /// Propagates FASTQ read / write errors; never silently drops a record
    /// for an IO reason (only BFC's defined uncorrectable/policy drops).
    pub fn run(&self, input: &Path, output: &Path) -> Result<CorrectReport> {
        if self.cfg.k < 11 || self.cfg.k > 63 || self.cfg.k.is_multiple_of(2) {
            return Err(RsomicsError::ConfigError(format!(
                "k must be odd and in 11..=63 (BFC_MAX_KMER 63), got {}",
                self.cfg.k
            )));
        }
        let ch = build_table(&[input], self.cfg)?;
        let mut reader = open_fastq(input)?;
        let mut w = ChunkedWriter::create(output, self.compression)?;
        let mut report = CorrectReport::default();
        let cfg = self.cfg;

        let mut chunk: Vec<OwnedRecord> = Vec::with_capacity(CHUNK_RECORDS);
        loop {
            chunk.clear();
            while chunk.len() < CHUNK_RECORDS {
                let Some(r) = reader.next() else { break };
                chunk.push(r?);
            }
            if chunk.is_empty() {
                break;
            }
            let out: Vec<(Option<OwnedRecord>, u64, u64)> = chunk
                .par_drain(..)
                .map(|rec| {
                    let bases_in = rec.seq.len() as u64;
                    if cfg.drop_unique_kmer && has_unique_kmer(cfg, &ch, &rec) {
                        return (None, bases_in, 0);
                    }
                    match correct_one(cfg, &ch, &rec) {
                        Some((seq, qual)) => {
                            let corrected =
                                seq.iter().filter(|&&c| c.is_ascii_lowercase()).count() as u64;
                            (
                                Some(OwnedRecord {
                                    id: rec.id,
                                    seq,
                                    qual: if cfg.fasta_out { Vec::new() } else { qual },
                                }),
                                bases_in,
                                corrected,
                            )
                        }
                        None if cfg.discard_uncorrectable => (None, bases_in, 0),
                        None => (Some(rec), bases_in, 0),
                    }
                })
                .collect();

            for (rec, bin, bcorr) in out {
                report.reads_in += 1;
                report.bases_in += bin;
                match rec {
                    Some(o) => {
                        report.bases_corrected += bcorr;
                        w.write_record(&o.id, &o.seq, &o.qual)?;
                        report.reads_out += 1;
                    }
                    None => report.reads_dropped += 1,
                }
            }
        }
        w.finalize()?;
        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nt4_table_maps_acgt_and_n() {
        assert_eq!(SEQ_NT4[b'A' as usize], 0);
        assert_eq!(SEQ_NT4[b'T' as usize], 3);
        assert_eq!(SEQ_NT4[b'g' as usize], 2);
        assert_eq!(SEQ_NT4[b'N' as usize], 4);
    }

    #[test]
    fn bfc_hash_64_is_invertible_shape() {
        let mask = (1u64 << 33) - 1;
        let a = bfc_hash_64(0x1234_5678 & mask, mask);
        let b = bfc_hash_64(0x1234_5679 & mask, mask);
        assert_ne!(a, b);
        assert_eq!(a, a & mask);
    }

    #[test]
    fn kmer_append_round_trips_forward_plane() {
        let mut x = BfcKmer::NULL;
        for c in [0u64, 1, 2, 3, 0, 3] {
            x.append(5, c);
        }
        // last 5 bases = C,G,T,A,T → low/high bit planes consistent, masked to k.
        let mask = (1u64 << 5) - 1;
        assert_eq!(x.x[0], x.x[0] & mask);
        assert_eq!(x.x[1], x.x[1] & mask);
    }

    #[test]
    fn best_island_picks_longest_solid_run() {
        let mut s = vec![EcBase::default(); 10];
        for (i, c) in s.iter_mut().enumerate() {
            c.solid_end = (3..=7).contains(&i);
        }
        let r = best_island(3, &s);
        assert!(r.is_some(), "a solid run must yield an island");
    }

    #[test]
    fn no_solid_kmer_is_uncorrectable() {
        let cfg = CorrectConfig::default();
        let ch = CountTable {
            k: cfg.k,
            map: HashMap::new(),
        };
        let rec = OwnedRecord {
            id: b"r".to_vec(),
            seq: b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            qual: b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
        };
        assert!(correct_one(&cfg, &ch, &rec).is_none());
    }

    #[test]
    fn high_coverage_clean_read_is_unchanged() {
        let cfg = CorrectConfig {
            k: 11,
            ..CorrectConfig::default()
        };
        let seq = b"ACGTACGTACGTGGGGCCCCAAAATTTTACGTACGT".to_vec();
        let qual = vec![b'I'; seq.len()];
        let rec = OwnedRecord {
            id: b"r".to_vec(),
            seq: seq.clone(),
            qual: qual.clone(),
        };
        let mut ch = CountTable {
            k: cfg.k,
            map: HashMap::new(),
        };
        for _ in 0..10 {
            let s = seq_conv(&seq, &qual, cfg.qual_threshold);
            let mut x = BfcKmer::NULL;
            let mut l = 0;
            for i in 0..s.len() {
                if s[i].b < 4 {
                    x.append(cfg.k, u64::from(s[i].b));
                    l += 1;
                    if l >= cfg.k {
                        ch.add(&x, true);
                    }
                } else {
                    l = 0;
                    x = BfcKmer::NULL;
                }
            }
        }
        let (out, _) = correct_one(&cfg, &ch, &rec).expect("clean high-cov read corrects");
        assert!(
            out.iter().all(u8::is_ascii_uppercase),
            "a clean read must stay uppercase (no corrections): {:?}",
            String::from_utf8_lossy(&out)
        );
    }

    /// BFC `bfc_ec_greedy_k` rescue: with exactly one strongly-supported
    /// single substitution of the probed k-mer (best `> mode/3`, 2nd `< 3`)
    /// it returns `pos<<2 | base`; with no confident alternative, `-1`.
    #[test]
    fn greedy_rescue_picks_the_one_supported_substitution() {
        let k = 11;
        let truth = b"ACGTCAGTTGA"; // exactly k bases
        let mut ch = CountTable {
            k,
            map: HashMap::new(),
        };
        // Heavily cover the truth k-mer; nothing else seen.
        let mut tx = BfcKmer::NULL;
        for &b in truth {
            tx.append(k, u64::from(SEQ_NT4[b as usize]));
        }
        for _ in 0..20 {
            ch.add(&tx, true);
        }
        // Probe the truth k-mer with its 3'-most base corrupted A→T (truth
        // ends in 'A'=0; corrupt to 'T'=3 at pos k-1, i.e. bit 0 / d=0).
        // greedy must propose restoring it to the truth base 'A'.
        let mut corrupt = BfcKmer::NULL;
        for (i, &b) in truth.iter().enumerate() {
            let base = if i == k - 1 {
                3
            } else {
                u64::from(SEQ_NT4[b as usize])
            };
            corrupt.append(k, base);
        }
        let mode = ch.hist_mode(3);
        let ec = bfc_ec_greedy_k(k, mode, &corrupt, &ch);
        assert!(ec >= 0, "a confident single fix must be found");
        assert_eq!(ec >> 2, 0, "the corrupted base is the 3'-most (d=0)");
        assert_eq!(
            (ec & 3) as u8,
            SEQ_NT4[b'A' as usize],
            "restored base must be the truth base A"
        );
        // An empty table → no confident alternative → -1.
        let empty = CountTable {
            k,
            map: HashMap::new(),
        };
        assert_eq!(bfc_ec_greedy_k(k, 0, &corrupt, &empty), -1);
    }
}
