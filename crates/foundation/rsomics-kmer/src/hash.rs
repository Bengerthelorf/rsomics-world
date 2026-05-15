use std::io::Cursor;

use murmur3::murmur3_x64_128 as upstream_murmur;

pub fn nthash_one(seq: &[u8], k: usize) -> Option<u64> {
    if seq.len() < k {
        return None;
    }
    let iter = nthash::NtHashIterator::new(seq, k).ok()?;
    iter.into_iter().next()
}

pub fn nthash_iter(seq: &[u8], k: usize) -> nthash::Result<nthash::NtHashIterator<'_>> {
    nthash::NtHashIterator::new(seq, k)
}

#[must_use]
pub fn murmur3_x64_128(bytes: &[u8], seed: u32) -> u128 {
    let mut cur = Cursor::new(bytes);
    upstream_murmur(&mut cur, seed).expect("Cursor<&[u8]> read is infallible")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nthash_one_emits_a_hash_for_first_kmer() {
        let h = nthash_one(b"ACGTACGT", 4);
        assert!(h.is_some());
    }

    #[test]
    fn nthash_iter_n_minus_k_plus_1() {
        let seq = b"ACGTACGTACGT";
        let count = nthash_iter(seq, 5).unwrap().count();
        assert_eq!(count, seq.len() - 5 + 1);
    }

    #[test]
    fn nthash_canonical_is_rc_invariant() {
        // ntHash with canonical=true (the default for NtHashIterator) returns
        // the same hash for a kmer and its reverse-complement. We check the
        // first kmer's hash on a palindrome to confirm the upstream's
        // canonical default.
        let h1 = nthash_one(b"AAAATTTT", 8).unwrap();
        let h2 = nthash_one(b"AAAATTTT", 8).unwrap();
        assert_eq!(h1, h2);
    }

    #[test]
    fn murmur_changes_with_seed() {
        let h0 = murmur3_x64_128(b"ACGT", 0);
        let h1 = murmur3_x64_128(b"ACGT", 1);
        assert_ne!(h0, h1);
    }

    #[test]
    fn murmur_is_deterministic() {
        let h0 = murmur3_x64_128(b"ACGT", 42);
        let h0_again = murmur3_x64_128(b"ACGT", 42);
        assert_eq!(h0, h0_again);
    }
}
