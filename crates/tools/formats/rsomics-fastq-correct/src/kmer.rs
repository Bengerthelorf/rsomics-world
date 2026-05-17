// Port of kmer.h (lh3/bfc, MIT).

// BFC bfc_kmer_t: x[0]/x[1] = 2-bit forward k-mer (low/high plane), x[2]/x[3] = RC; k ≤ 63.
#[derive(Clone, Copy)]
pub(crate) struct BfcKmer {
    pub(crate) x: [u64; 4],
}

impl BfcKmer {
    pub(crate) const NULL: BfcKmer = BfcKmer { x: [0, 0, 0, 0] };

    #[inline]
    pub(crate) fn append(&mut self, k: usize, c: u64) {
        // k validated 11..=63 at boundary; k≥64 unreachable — no branch in hot path.
        let mask = (1u64 << k) - 1;
        self.x[0] = ((self.x[0] << 1) | (c & 1)) & mask;
        self.x[1] = ((self.x[1] << 1) | (c >> 1)) & mask;
        self.x[2] = (self.x[2] >> 1) | ((1u64 ^ (c & 1)) << (k - 1));
        self.x[3] = (self.x[3] >> 1) | ((1u64 ^ (c >> 1)) << (k - 1));
    }

    // BFC bfc_kmer_change: d counted from 3' end (0 ≤ d < k); updates forward + RC planes.
    #[inline]
    pub(crate) fn change(&mut self, k: usize, d: usize, c: u64) {
        let t = !(1u64 << d);
        self.x[0] = ((c & 1) << d) | (self.x[0] & t);
        self.x[1] = ((c >> 1) << d) | (self.x[1] & t);
        let t = !(1u64 << (k - 1 - d));
        self.x[2] = ((1u64 ^ (c & 1)) << (k - 1 - d)) | (self.x[2] & t);
        self.x[3] = ((1u64 ^ (c >> 1)) << (k - 1 - d)) | (self.x[3] & t);
    }
}

// BFC bfc_hash_64 — Thomas Wang's invertible integer hash, masked.
#[inline]
pub(crate) fn bfc_hash_64(mut key: u64, mask: u64) -> u64 {
    key = (!key).wrapping_add(key << 21) & mask;
    key ^= key >> 24;
    key = key.wrapping_add(key << 3).wrapping_add(key << 8) & mask;
    key ^= key >> 14;
    key = key.wrapping_add(key << 2).wrapping_add(key << 4) & mask;
    key ^= key >> 28;
    key = key.wrapping_add(key << 31) & mask;
    key
}

// BFC bfc_kmer_hash: middle base selects strand; reproduced verbatim for compat
// (compat gap: our plain map vs BFC's open-address bfc_ch has different collision profile).
#[inline]
pub(crate) fn bfc_kmer_hash(k: usize, x: &[u64; 4]) -> u64 {
    // k < 64 required for `<< k`; validated at boundary, assert rather than silent zero.
    debug_assert!(k < 64, "k must be < 64 (BFC_MAX_KMER); got {k}");
    let t = k >> 1;
    let u = usize::from((x[1] >> t & 1) > (x[3] >> t & 1));
    let mask = (1u64 << k) - 1;
    let h0 = bfc_hash_64((x[u << 1].wrapping_add(x[(u << 1) | 1])) & mask, mask);
    let h1 = bfc_hash_64(h0 ^ x[(u << 1) | 1], mask);
    ((h0 ^ h1) << k) | ((h0.wrapping_add(h1)) & mask)
}
