#![allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum FmError {
    #[error(
        "input contains the sentinel byte 0x00 at offset {pos} — strip or remap before indexing"
    )]
    SentinelInInput { pos: usize },
    #[error("empty input")]
    Empty,
}

pub type Result<T> = std::result::Result<T, FmError>;

const SENTINEL: u8 = 0;
const OCC_STRIDE: usize = 32;

#[derive(Debug, Clone)]
pub struct FmIndex {
    pub bwt: Vec<u8>,
    pub sa: Vec<i32>,
    c: [usize; 256],      // C[c] = count of bytes < c in bwt
    occ: Vec<[u32; 256]>, // sparse rank table, stride = OCC_STRIDE
}

impl FmIndex {
    pub fn build(text: &[u8]) -> Result<Self> {
        if text.is_empty() {
            return Err(FmError::Empty);
        }
        if let Some(pos) = text.iter().position(|&b| b == SENTINEL) {
            return Err(FmError::SentinelInInput { pos });
        }
        let mut t = Vec::with_capacity(text.len() + 1);
        t.extend_from_slice(text);
        t.push(SENTINEL);

        let mut sa = vec![0_i32; t.len()];
        divsufsort::sort_in_place(&t, &mut sa);

        let mut bwt = vec![0_u8; t.len()];
        for i in 0..t.len() {
            let pos = sa[i] as usize;
            bwt[i] = if pos == 0 { t[t.len() - 1] } else { t[pos - 1] };
        }

        let mut counts = [0_usize; 256];
        for &b in &bwt {
            counts[b as usize] += 1;
        }
        let mut c = [0_usize; 256];
        let mut acc = 0_usize;
        for (i, item) in c.iter_mut().enumerate() {
            *item = acc;
            acc += counts[i];
        }

        let n_blocks = bwt.len().div_ceil(OCC_STRIDE) + 1;
        let mut occ = vec![[0_u32; 256]; n_blocks];
        let mut running = [0_u32; 256];
        for i in 0..bwt.len() {
            if i % OCC_STRIDE == 0 {
                occ[i / OCC_STRIDE] = running;
            }
            running[bwt[i] as usize] += 1;
        }
        occ[n_blocks - 1] = running;

        Ok(Self { bwt, sa, c, occ })
    }

    #[must_use]
    pub fn occ(&self, b: u8, i: usize) -> usize {
        let block = i / OCC_STRIDE;
        let base = self.occ[block][b as usize] as usize;
        let mut acc = base;
        let block_start = block * OCC_STRIDE;
        for &x in &self.bwt[block_start..i] {
            if x == b {
                acc += 1;
            }
        }
        acc
    }

    /// SA interval `[lo, hi)` of suffixes starting with `pattern`, or `None`.
    #[must_use]
    pub fn backward_search(&self, pattern: &[u8]) -> Option<(usize, usize)> {
        if pattern.is_empty() {
            return Some((0, self.bwt.len()));
        }
        let mut lo = 0_usize;
        let mut hi = self.bwt.len();
        for &c in pattern.iter().rev() {
            let cc = self.c[c as usize];
            lo = cc + self.occ(c, lo);
            hi = cc + self.occ(c, hi);
            if lo >= hi {
                return None;
            }
        }
        Some((lo, hi))
    }

    #[must_use]
    pub fn count(&self, pattern: &[u8]) -> usize {
        self.backward_search(pattern).map_or(0, |(lo, hi)| hi - lo)
    }

    #[must_use]
    pub fn locate(&self, pattern: &[u8]) -> Vec<usize> {
        match self.backward_search(pattern) {
            Some((lo, hi)) => self.sa[lo..hi].iter().map(|&p| p as usize).collect(),
            None => Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn count_matches_in_simple_text() {
        let fm = FmIndex::build(b"abracadabra").unwrap();
        assert_eq!(fm.count(b"a"), 5);
        assert_eq!(fm.count(b"abra"), 2);
        assert_eq!(fm.count(b"rac"), 1);
        assert_eq!(fm.count(b"xyz"), 0);
    }

    #[test]
    fn locate_returns_all_positions() {
        let fm = FmIndex::build(b"abracadabra").unwrap();
        let mut positions = fm.locate(b"abra");
        positions.sort_unstable();
        assert_eq!(positions, vec![0, 7]);
    }

    #[test]
    fn locate_single_char_matches() {
        let fm = FmIndex::build(b"abracadabra").unwrap();
        let mut p = fm.locate(b"r");
        p.sort_unstable();
        assert_eq!(p, vec![2, 9]);
    }

    #[test]
    fn count_zero_on_missing_pattern() {
        let fm = FmIndex::build(b"acgtacgt").unwrap();
        assert_eq!(fm.count(b"ggg"), 0);
        assert_eq!(fm.locate(b"ggg"), Vec::<usize>::new());
    }

    #[test]
    fn empty_pattern_returns_full_sa() {
        let fm = FmIndex::build(b"acgt").unwrap();
        let (lo, hi) = fm.backward_search(b"").unwrap();
        assert_eq!((lo, hi), (0, fm.bwt.len()));
    }

    #[test]
    fn sentinel_in_input_rejected() {
        let mut text = b"acgt".to_vec();
        text.push(0);
        text.extend_from_slice(b"acgt");
        assert!(matches!(
            FmIndex::build(&text),
            Err(FmError::SentinelInInput { pos: 4 })
        ));
    }

    #[test]
    fn empty_input_rejected() {
        assert!(matches!(FmIndex::build(b""), Err(FmError::Empty)));
    }

    #[test]
    fn dna_pattern_search() {
        let fm = FmIndex::build(b"AAGATCAATGATCAAAA").unwrap();
        assert_eq!(fm.count(b"GATC"), 2);
        let mut p = fm.locate(b"GATC");
        p.sort_unstable();
        assert_eq!(p, vec![2, 9]);
    }
}
