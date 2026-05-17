use crate::encode::{Kmer, base_bits};

pub struct RollingKmers<'a> {
    seq: &'a [u8],
    k: usize,
    mask: u64,
    pos: usize,
    current: u64,
    valid: usize,
}

impl<'a> RollingKmers<'a> {
    #[must_use]
    pub fn new(seq: &'a [u8], k: usize) -> Self {
        debug_assert!((1..=32).contains(&k));
        Self {
            seq,
            k,
            mask: if k == 32 {
                u64::MAX
            } else {
                (1u64 << (2 * k)) - 1
            },
            pos: 0,
            current: 0,
            valid: 0,
        }
    }
}

impl Iterator for RollingKmers<'_> {
    type Item = Option<Kmer>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.seq.len() {
            return None;
        }
        let b = self.seq[self.pos];
        self.pos += 1;

        if let Some(bits) = base_bits(b) {
            self.current = ((self.current << 2) | bits) & self.mask;
            self.valid += 1;
        } else {
            self.current = 0;
            self.valid = 0;
        }

        if self.valid >= self.k {
            Some(Some(self.current))
        } else {
            Some(None)
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.seq.len() - self.pos;
        (remaining, Some(remaining))
    }
}

impl ExactSizeIterator for RollingKmers<'_> {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode::encode;

    #[test]
    fn rolling_matches_encode() {
        let seq = b"ACGTACGTACGT";
        let k = 4;
        let rolling: Vec<u64> = RollingKmers::new(seq, k).flatten().collect();
        let naive: Vec<u64> = seq.windows(k).map(|w| encode(w).unwrap()).collect();
        assert_eq!(rolling, naive);
    }

    #[test]
    fn rolling_skips_n_bearing_windows() {
        let seq = b"ACGTNACGT";
        let k = 4;
        let results: Vec<Option<u64>> = RollingKmers::new(seq, k).collect();
        assert_eq!(results.len(), 9);
        // First valid k-mer at index k-1=3 (ACGT)
        assert!(results[0].is_none());
        assert!(results[1].is_none());
        assert!(results[2].is_none());
        assert!(results[3].is_some()); // ACGT
        // N at index 4 resets; need 4 more valid bases to recover
        assert!(results[4].is_none()); // N resets
        assert!(results[5].is_none()); // A (valid=1)
        assert!(results[6].is_none()); // C (valid=2)
        assert!(results[7].is_none()); // G (valid=3)
        assert!(results[8].is_some()); // T (valid=4 → ACGT)
    }

    #[test]
    fn rolling_empty_seq() {
        let results: Vec<_> = RollingKmers::new(b"", 4).collect();
        assert!(results.is_empty());
    }

    #[test]
    fn rolling_seq_shorter_than_k() {
        let results: Vec<_> = RollingKmers::new(b"ACG", 4).collect();
        assert_eq!(results.len(), 3);
        assert!(results.iter().all(Option::is_none));
    }
}
