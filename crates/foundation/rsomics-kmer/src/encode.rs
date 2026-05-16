use crate::{KmerError, Result};

pub type Kmer = u64;

const MAX_K: usize = 32;

const fn base_bits(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0b00),
        b'C' | b'c' => Some(0b01),
        b'G' | b'g' => Some(0b10),
        b'T' | b't' => Some(0b11),
        _ => None,
    }
}

const fn bits_base(bits: u64) -> u8 {
    match bits & 0b11 {
        0b00 => b'A',
        0b01 => b'C',
        0b10 => b'G',
        _ => b'T',
    }
}

pub fn encode(seq: &[u8]) -> Result<Kmer> {
    let k = seq.len();
    if !(1..=MAX_K).contains(&k) {
        return Err(KmerError::KOutOfRange(k));
    }
    let mut bits: u64 = 0;
    for (i, &b) in seq.iter().enumerate() {
        let v = base_bits(b).ok_or(KmerError::NonAcgt { pos: i, byte: b })?;
        bits = (bits << 2) | v;
    }
    Ok(bits)
}

#[must_use]
pub fn decode(kmer: Kmer, k: usize) -> Vec<u8> {
    let mut out = vec![0u8; k];
    let mut bits = kmer;
    for slot in out.iter_mut().rev() {
        *slot = bits_base(bits);
        bits >>= 2;
    }
    out
}

#[must_use]
pub fn reverse_complement(kmer: Kmer, k: usize) -> Kmer {
    let mut bits = kmer;
    let comp = bits ^ ((1u64 << (2 * k)) - 1); // XOR with all-1s per 2-bit pair: A↔T, C↔G
    let mut rc: u64 = 0;
    bits = comp;
    for _ in 0..k {
        rc = (rc << 2) | (bits & 0b11);
        bits >>= 2;
    }
    rc
}

#[must_use]
pub fn canonical(kmer: Kmer, k: usize) -> Kmer {
    let rc = reverse_complement(kmer, k);
    kmer.min(rc)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_round_trip() {
        let k = 8;
        let seq = b"ACGTACGT";
        let bits = encode(seq).unwrap();
        let back = decode(bits, k);
        assert_eq!(&back, seq);
    }

    #[test]
    fn encode_rejects_n() {
        let r = encode(b"ACGTN");
        assert!(matches!(r, Err(KmerError::NonAcgt { pos: 4, .. })));
    }

    #[test]
    fn encode_k_too_large_rejected() {
        let seq = vec![b'A'; 33];
        assert!(matches!(encode(&seq), Err(KmerError::KOutOfRange(33))));
    }

    #[test]
    fn rc_of_rc_is_identity() {
        let bits = encode(b"ACGTACGT").unwrap();
        let rc = reverse_complement(bits, 8);
        let back = reverse_complement(rc, 8);
        assert_eq!(back, bits);
    }

    #[test]
    fn rc_known_value() {
        assert_eq!(
            decode(reverse_complement(encode(b"AAAA").unwrap(), 4), 4),
            b"TTTT".to_vec()
        );
        assert_eq!(
            decode(reverse_complement(encode(b"ACGT").unwrap(), 4), 4),
            b"ACGT".to_vec()
        );
    }

    #[test]
    fn canonical_picks_lex_min_of_pair() {
        let fwd = encode(b"GGGG").unwrap();
        let rc = reverse_complement(fwd, 4);
        assert_eq!(canonical(fwd, 4), rc);
        assert_eq!(canonical(fwd, 4), canonical(rc, 4));
    }
}
