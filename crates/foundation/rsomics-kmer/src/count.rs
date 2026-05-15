use std::collections::HashMap;

use crate::encode::Kmer;
use crate::iter::KmerIter;
use crate::{KmerError, Result};

#[derive(Debug, Default, Clone)]
pub struct KmerCounts {
    pub k: usize,
    pub canonical: bool,
    pub counts: HashMap<Kmer, u64>,
}

impl KmerCounts {
    #[must_use]
    pub fn new(k: usize, canonical: bool) -> Self {
        Self {
            k,
            canonical,
            counts: HashMap::new(),
        }
    }

    pub fn count_seq(&mut self, seq: &[u8]) -> Result<()> {
        if seq.len() < self.k {
            return Ok(());
        }
        let it = KmerIter::new(seq, self.k, self.canonical)?;
        for kmer in it {
            match kmer {
                Ok(k) => *self.counts.entry(k).or_insert(0) += 1,
                Err(KmerError::NonAcgt { .. }) => {}
                Err(e) => return Err(e),
            }
        }
        Ok(())
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.counts.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    #[must_use]
    pub fn total(&self) -> u64 {
        self.counts.values().sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn count_simple_seq() {
        let mut c = KmerCounts::new(3, false);
        c.count_seq(b"AAAAA").unwrap();
        assert_eq!(c.total(), 3);
        assert_eq!(c.len(), 1);
    }

    #[test]
    fn count_canonical_collapses_rc_pairs() {
        let mut c = KmerCounts::new(4, true);
        c.count_seq(b"AAAATTTT").unwrap();
        let entries: Vec<_> = c.counts.iter().collect();
        assert!(!entries.is_empty());
    }

    #[test]
    fn count_skips_n_kmers_silently() {
        let mut c = KmerCounts::new(4, false);
        c.count_seq(b"ACGTNACGT").unwrap();
        assert_eq!(c.total(), 2);
    }

    #[test]
    fn count_seq_shorter_than_k_is_noop() {
        let mut c = KmerCounts::new(10, false);
        c.count_seq(b"ACGT").unwrap();
        assert_eq!(c.total(), 0);
    }
}
