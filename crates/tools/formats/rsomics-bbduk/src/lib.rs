// clap --help text names BBDuk / bbduk.sh params that read as code; backticking them clutters the help
#![allow(clippy::doc_markdown)]

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

use rsomics_kmer::{RollingKmers, encode, reverse_complement};

// k-mer codes are already well-distributed 2-bit encodings; SipHash is wasted work
#[derive(Default)]
struct IdentityHasher(u64);

impl Hasher for IdentityHasher {
    fn write_u64(&mut self, n: u64) {
        self.0 = n;
    }
    fn write(&mut self, _: &[u8]) {
        unreachable!();
    }
    fn finish(&self) -> u64 {
        self.0
    }
}

type KmerSet = std::collections::HashSet<u64, BuildHasherDefault<IdentityHasher>>;

// k ≤ 31: reverse_complement does 1u64 << (2*k), UB at k=32; BBDuk's own default is 27
pub const MAX_K: usize = 31;

// variant index grows (3k)^hdist; past 3 it outgrows memory — fail loud, not OOM
pub const MAX_HDIST: usize = 3;

// BBDuk ktrim: None=f (kfilter, a hit removes the read), Right=r (trim hit→3'), Left=l (trim 5'→hit end)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KTrim {
    None,
    Right,
    Left,
}

// BBDuk qtrim=f|r|l|rl
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QTrim {
    None,
    Right,
    Left,
    Both,
}

// fields + defaults track bbduk.sh's documented defaults (see impl Default)
#[derive(Debug, Clone)]
pub struct Config {
    pub k: usize,
    // BBDuk mink: also index/probe ref k-mers this short at read tips (0 = off) for partial-adapter trim
    pub mink: usize,
    // BBDuk hdist: index every ref k-mer + all k-mers within this Hamming distance (0 = exact)
    pub hdist: usize,
    // BBDuk hdist2: Hamming distance for the short (mink) tip k-mers — kept separate from hdist and
    // defaulting to 0, so --hdist 1 --mink 11 means fuzzy full k-mers but exact tips (faithful, not a bug)
    pub hdist2: usize,
    // BBDuk rcomp: also match each ref k-mer's reverse complement
    pub rcomp: bool,
    // BBDuk maskmiddle: wildcard each k-mer's middle base; auto-disabled when mink>0 (RefKmers::build enforces)
    pub maskmiddle: bool,
    // BBDuk minkmerhits (mkh): contaminant if it shares ≥ this many k-mers with the reference
    pub min_kmer_hits: usize,
    // BBDuk minkmerfraction (mkf): …or ≥ this fraction of its k-mers
    pub min_kmer_fraction: f64,
    pub ktrim: KTrim,
    pub qtrim: QTrim,
    // BBDuk trimq: quality-trim threshold (Phred)
    pub trimq: u8,
    // BBDuk minlength: discard reads shorter than this after trimming
    pub min_length: usize,
    pub qual_offset: u8,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            k: 27,
            mink: 0,
            hdist: 0,
            hdist2: 0,
            rcomp: true,
            maskmiddle: true,
            min_kmer_hits: 1,
            min_kmer_fraction: 0.0,
            ktrim: KTrim::None,
            qtrim: QTrim::None,
            trimq: 6,
            min_length: 10,
            qual_offset: 33,
        }
    }
}

// zero the middle base's 2-bit field (BBDuk maskmiddle); for odd k the centre is revcomp-invariant, so masking at insert and probe is symmetric
#[must_use]
fn mask_mid(code: u64, k: usize) -> u64 {
    let m = k / 2;
    let shift = 2 * (k - 1 - m);
    code & !(0b11u64 << shift)
}

// non-ACGT window → nothing: BBDuk drops ambiguous ref k-mers; expanding one would fabricate reference content
fn variants(w: &[u8], hdist: usize) -> Vec<u64> {
    if w.iter().any(|b| !matches!(b, b'A' | b'C' | b'G' | b'T')) {
        return Vec::new();
    }
    let mut buf = w.to_vec();
    let mut out = Vec::new();
    expand(&mut buf, 0, hdist, &mut out);
    out
}

fn expand(buf: &mut Vec<u8>, start: usize, budget: usize, out: &mut Vec<u64>) {
    if let Ok(code) = encode(buf) {
        out.push(code);
    }
    if budget == 0 {
        return;
    }
    for pos in start..buf.len() {
        let orig = buf[pos];
        for &b in b"ACGT" {
            if b != orig {
                buf[pos] = b;
                expand(buf, pos + 1, budget - 1, out);
            }
        }
        buf[pos] = orig;
    }
}

// exact ref k-mer index (not a Bloom filter — false positives would break kfilter byte-compat);
// mink>0 also indexes shorter tip k-mers (lengths mink..k) governed by hdist2
pub struct RefKmers {
    full: KmerSet,
    tips: HashMap<usize, KmerSet>, // length → ref k-mers of that length (mink..k)
    k: usize,
    mink: usize,
    // effective maskmiddle = requested && mink==0 — BBDuk disables maskmiddle when short k-mers are in play
    mm: bool,
}

impl RefKmers {
    #[must_use]
    pub fn build<'a>(refs: impl IntoIterator<Item = &'a [u8]>, cfg: &Config) -> Self {
        assert!(
            (1..=MAX_K).contains(&cfg.k),
            "k must be in 1..={MAX_K} (2-bit codec / rcomp limit); got {}",
            cfg.k
        );
        assert!(
            cfg.hdist <= MAX_HDIST && cfg.hdist2 <= MAX_HDIST,
            "hdist/hdist2 must be ≤ {MAX_HDIST} (variant index memory bound)"
        );
        assert!(cfg.mink <= cfg.k, "mink must be ≤ k");

        let mm = cfg.maskmiddle && cfg.mink == 0;
        let key = |c: u64| if mm { mask_mid(c, cfg.k) } else { c };

        let mut full = KmerSet::default();
        let mut tips: HashMap<usize, KmerSet> = HashMap::new();
        let tip_lens: Vec<usize> = if cfg.mink > 0 {
            (cfg.mink..cfg.k).collect()
        } else {
            Vec::new()
        };
        for seq in refs {
            for w in seq.windows(cfg.k) {
                for code in variants(w, cfg.hdist) {
                    full.insert(key(code));
                    if cfg.rcomp {
                        full.insert(key(reverse_complement(code, cfg.k)));
                    }
                }
            }
            for &l in &tip_lens {
                if seq.len() < l {
                    continue;
                }
                let set = tips.entry(l).or_default();
                for w in seq.windows(l) {
                    for code in variants(w, cfg.hdist2) {
                        set.insert(code);
                        if cfg.rcomp {
                            set.insert(reverse_complement(code, l));
                        }
                    }
                }
            }
        }
        Self {
            full,
            tips,
            k: cfg.k,
            mink: cfg.mink,
            mm,
        }
    }

    fn key(&self, c: u64) -> u64 {
        if self.mm { mask_mid(c, self.k) } else { c }
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.full.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.full.is_empty()
    }

    #[must_use]
    pub fn hits(&self, seq: &[u8]) -> usize {
        if seq.len() < self.k {
            return 0;
        }
        RollingKmers::new(seq, self.k)
            .flatten()
            .filter(|&c| self.full.contains(&self.key(c)))
            .count()
    }

    #[must_use]
    pub fn first_full_hit(&self, seq: &[u8]) -> Option<usize> {
        if seq.len() < self.k {
            return None;
        }
        RollingKmers::new(seq, self.k)
            .enumerate()
            .find_map(|(i, opt)| {
                opt.filter(|&c| self.full.contains(&self.key(c)))
                    .map(|_| i + 1 - self.k)
            })
    }

    // ktrim=r cut: full-k-mer hit, else (mink>0) start of the longest ref-matching 3' tip — BBDuk never trims an internal short k-mer
    #[must_use]
    pub fn right_cut(&self, seq: &[u8]) -> Option<usize> {
        if let Some(p) = self.first_full_hit(seq) {
            return Some(p);
        }
        if self.mink == 0 {
            return None;
        }
        // when the read is shorter than k it is itself the longest tip, so the bound is inclusive: min(k-1, seq.len())
        for l in (self.mink..=(self.k - 1).min(seq.len())).rev() {
            let tip = &seq[seq.len() - l..];
            if self
                .tips
                .get(&l)
                .is_some_and(|s| encode(tip).is_ok_and(|c| s.contains(&c)))
            {
                return Some(seq.len() - l);
            }
        }
        None
    }

    // ktrim=l keep-from: end of full-k-mer hit, else (mink>0) end of the longest ref-matching 5' tip
    #[must_use]
    pub fn left_cut(&self, seq: &[u8]) -> Option<usize> {
        if let Some(p) = self.first_full_hit(seq) {
            return Some(p + self.k);
        }
        if self.mink == 0 {
            return None;
        }
        for l in (self.mink..=(self.k - 1).min(seq.len())).rev() {
            let tip = &seq[..l];
            if self
                .tips
                .get(&l)
                .is_some_and(|s| encode(tip).is_ok_and(|c| s.contains(&c)))
            {
                return Some(l);
            }
        }
        None
    }
}

// BBDuk kfilter contaminant test
#[must_use]
pub fn is_contaminant(seq: &[u8], refs: &RefKmers, cfg: &Config) -> bool {
    let hits = refs.hits(seq);
    if hits == 0 {
        return false;
    }
    if hits >= cfg.min_kmer_hits {
        return true;
    }
    if cfg.min_kmer_fraction > 0.0 {
        let total = seq.len().saturating_sub(cfg.k) + 1;
        if total > 0 {
            #[allow(clippy::cast_precision_loss)]
            let frac = hits as f64 / total as f64;
            return frac >= cfg.min_kmer_fraction;
        }
    }
    false
}

// BBDuk pipeline order: k-mer filter / k-mer trim → quality trim → min-length
#[must_use]
pub fn process(seq: &[u8], qual: &[u8], refs: &RefKmers, cfg: &Config) -> Option<(usize, usize)> {
    let mut start = 0usize;
    let mut end = seq.len();

    if !refs.is_empty() {
        match cfg.ktrim {
            KTrim::None => {
                if is_contaminant(seq, refs, cfg) {
                    return None;
                }
            }
            KTrim::Right => {
                if let Some(pos) = refs.right_cut(seq) {
                    end = pos;
                }
            }
            KTrim::Left => {
                if let Some(pos) = refs.left_cut(seq) {
                    start = pos.min(end);
                }
            }
        }
    }

    if end > start && cfg.qtrim != QTrim::None {
        let thr = cfg.trimq;
        let phred = |q: u8| q.saturating_sub(cfg.qual_offset);
        if matches!(cfg.qtrim, QTrim::Right | QTrim::Both) {
            while end > start && phred(qual[end - 1]) < thr {
                end -= 1;
            }
        }
        if matches!(cfg.qtrim, QTrim::Left | QTrim::Both) {
            while start < end && phred(qual[start]) < thr {
                start += 1;
            }
        }
    }

    if end.saturating_sub(start) < cfg.min_length {
        return None;
    }
    Some((start, end))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn refs(seqs: &[&[u8]], cfg: &Config) -> RefKmers {
        RefKmers::build(seqs.iter().copied(), cfg)
    }

    // Exact-match property tests disable maskmiddle so the centre base
    // is not wildcarded.
    fn exact(k: usize) -> Config {
        Config {
            k,
            maskmiddle: false,
            ..Config::default()
        }
    }

    #[test]
    fn contaminant_read_with_shared_kmer_is_filtered() {
        let cfg = exact(11);
        let r = refs(&[b"ACGTTGCAACGTTGCA"], &cfg);
        assert!(is_contaminant(b"TTTTACGTTGCAACGTTGCATTTT", &r, &cfg));
        assert!(!is_contaminant(b"TTTTTTTTTTTTTTTTTTTTTTTT", &r, &cfg));
    }

    #[test]
    fn rcomp_match_is_found_when_enabled() {
        let fwd: &[u8] = b"ACGTTGCAACG";
        let on = exact(11);
        let off = Config {
            rcomp: false,
            ..on.clone()
        };
        let probe: Vec<u8> = fwd
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                _ => b'C',
            })
            .collect();
        assert!(is_contaminant(&probe, &refs(&[fwd], &on), &on));
        assert!(!is_contaminant(&probe, &refs(&[fwd], &off), &off));
    }

    #[test]
    fn hdist_one_matches_single_mismatch_kmer() {
        let cfg = Config {
            hdist: 1,
            rcomp: false,
            ..exact(11)
        };
        let r = refs(&[b"ACGTTGCAACG"], &cfg);
        // index 1 is not the middle (k=11 centre is index 5), so hdist=0 must reject it.
        assert!(is_contaminant(b"ATGTTGCAACG", &r, &cfg));
        let strict = Config {
            hdist: 0,
            ..cfg.clone()
        };
        assert!(!is_contaminant(
            b"ATGTTGCAACG",
            &refs(&[b"ACGTTGCAACG"], &strict),
            &strict
        ));
    }

    #[test]
    fn maskmiddle_default_wildcards_the_centre_base() {
        // k=11 ⇒ middle index 5. A single mismatch *at the middle* is a
        // hit with maskmiddle (the BBDuk default), a miss without it.
        let mm = Config {
            k: 11,
            ..Config::default()
        };
        let r = refs(&[b"ACGTTGCAACG"], &mm);
        assert!(is_contaminant(b"ACGTTACAACG", &r, &mm)); // index-5 G→A
        let no_mm = exact(11);
        assert!(!is_contaminant(
            b"ACGTTACAACG",
            &refs(&[b"ACGTTGCAACG"], &no_mm),
            &no_mm
        ));
    }

    #[test]
    fn mink_disables_maskmiddle() {
        // mink>0 overrides the defaulted-true maskmiddle.
        let cfg = Config {
            k: 11,
            mink: 6,
            ..Config::default()
        };
        let r = refs(&[b"ACGTTGCAACG"], &cfg);
        assert!(!is_contaminant(b"ACGTTACAACG", &r, &cfg));
    }

    #[test]
    fn ktrim_right_trims_from_adapter_kmer() {
        let cfg = Config {
            ktrim: KTrim::Right,
            rcomp: false,
            ..exact(12)
        };
        let adapter = b"AGATCGGAAGAGC";
        let r = refs(&[adapter], &cfg);
        let read = b"ACGTTGCAACGTACGTTGCAAGATCGGAAGAGC";
        let qual = vec![b'I'; read.len()];
        let (s, e) = process(read, &qual, &r, &cfg).unwrap();
        assert_eq!(s, 0);
        assert_eq!(e, 20);
    }

    #[test]
    fn mink_trims_partial_adapter_at_read_tip() {
        let cfg = Config {
            k: 12,
            mink: 6,
            ktrim: KTrim::Right,
            rcomp: false,
            ..Config::default()
        };
        let adapter = b"AGATCGGAAGAGCACAC";
        let r = refs(&[adapter], &cfg);
        // only the first 7 bp of the adapter are present (< k=12, ≥ mink=6).
        let read = b"ACGTTGCAACGTACGTTGCAAGATCGG";
        let qual = vec![b'I'; read.len()];
        let (s, e) = process(read, &qual, &r, &cfg).unwrap();
        assert_eq!(s, 0);
        assert_eq!(e, 20);
    }

    #[test]
    fn short_read_whole_length_tip_is_checked() {
        // The tip-length upper bound is inclusive: a read shorter than k
        // whose entire body is the adapter tip must still clip.
        let cfg = Config {
            k: 12,
            mink: 6,
            ktrim: KTrim::Right,
            rcomp: false,
            ..Config::default()
        };
        let adapter = b"AGATCGGAAGAGCACAC";
        let r = refs(&[adapter], &cfg);
        let read = b"AGATCGG";
        let qual = vec![b'I'; read.len()];
        // entirely adapter ⇒ trimmed to length 0 ⇒ discarded (< minlength).
        assert!(process(read, &qual, &r, &cfg).is_none());
    }

    #[test]
    fn variants_skip_ambiguous_ref_window() {
        // An N in the ref window must not fabricate ACGT variants.
        let cfg = Config {
            hdist: 1,
            rcomp: false,
            ..exact(11)
        };
        let r = refs(&[b"ACGTNGCAACG"], &cfg);
        assert!(r.is_empty(), "N-bearing ref window must index nothing");
    }

    #[test]
    fn qtrim_right_and_minlength_discard() {
        let cfg = Config {
            qtrim: QTrim::Right,
            trimq: 20,
            min_length: 10,
            ..Config::default()
        };
        let empty = RefKmers::build(std::iter::empty(), &cfg);
        let seq = b"ACGTACGTACGTACGTACGT";
        let mut q = vec![b'I'; 8];
        q.extend(std::iter::repeat_n(b'#', 12));
        assert!(process(seq, &q, &empty, &cfg).is_none());
    }

    #[test]
    fn clean_read_survives_unchanged() {
        let cfg = Config::default();
        let empty = RefKmers::build(std::iter::empty(), &cfg);
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let q = vec![b'I'; seq.len()];
        assert_eq!(process(seq, &q, &empty, &cfg), Some((0, seq.len())));
    }
}
