// Port of htab.c/bbf.c (lh3/bfc, MIT).

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::path::Path;

use rsomics_common::Result;
use rsomics_seqio::{OwnedRecord, open_fastq};

use crate::CorrectConfig;
use crate::correct::seq_conv;
use crate::kmer::{BfcKmer, bfc_kmer_hash};

// Key is already a Thomas-Wang double-hash (bfc_kmer_hash); SipHash re-hash was ~half the
// hot path. BFC's bfc_ch also indexes directly by hash. Only u64 keys are valid.
#[derive(Default)]
pub(crate) struct IdentityHasher(u64);

impl Hasher for IdentityHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, _: &[u8]) {
        unreachable!("CountTable key is u64; only write_u64 is used");
    }
    fn write_u64(&mut self, n: u64) {
        self.0 = n;
    }
}

pub(crate) type KmerMap = HashMap<u64, Occ, BuildHasherDefault<IdentityHasher>>;

// BFC bfc_ch_kmer_occ layout: cov = r&0xff, hi = r>>8&0x3f; both saturating.
#[derive(Default, Clone, Copy)]
pub(crate) struct Occ {
    pub(crate) cov: u8,
    pub(crate) hi: u8,
}

pub(crate) struct CountTable {
    pub(crate) k: usize,
    pub(crate) map: KmerMap,
}

impl CountTable {
    pub(crate) fn occ(&self, kmer: &BfcKmer) -> Occ {
        self.map
            .get(&bfc_kmer_hash(self.k, &kmer.x))
            .copied()
            .unwrap_or_default()
    }

    pub(crate) fn add(&mut self, kmer: &BfcKmer, high_qual: bool) {
        let e = self.map.entry(bfc_kmer_hash(self.k, &kmer.x)).or_default();
        e.cov = e.cov.saturating_add(1);
        if high_qual {
            e.hi = e.hi.saturating_add(1).min(0x3f);
        }
    }

    // BFC bfc_ch_hist: peak over cov ≥ min_cov to skip the error-noise spike at low coverage.
    pub(crate) fn hist_mode(&self, min_cov: i32) -> i32 {
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

pub(crate) fn build_table(paths: &[&Path], cfg: &CorrectConfig) -> Result<CountTable> {
    let mut ch = CountTable {
        k: cfg.k,
        map: KmerMap::default(),
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

pub(crate) fn has_unique_kmer(cfg: &CorrectConfig, ch: &CountTable, rec: &OwnedRecord) -> bool {
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
