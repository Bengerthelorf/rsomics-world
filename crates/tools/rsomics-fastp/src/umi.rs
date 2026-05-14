//! UMI (Unique Molecular Identifier) extraction.
//!
//! Takes the first `len` bases of a read as the UMI, appends them to the
//! read id (fastp's convention: `:` separator), and returns the trimmed
//! seq + qual ready for the downstream filter pipeline. The subset we expose
//! is `umi_loc=read1` / `umi_loc=read2` (UMI embedded in a sequenced mate);
//! the `per_index` case (UMI in a separate index file) is not part of this
//! module's scope.

/// Where to pull the UMI bases from. `Read1` extracts from R1; `Read2` from
/// R2. For SE only `Read1` is valid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmiLoc {
    Read1,
    Read2,
}

#[derive(Debug, Clone, Copy)]
pub struct UmiConfig {
    pub loc: UmiLoc,
    pub len: usize,
}

/// Extract the first `cfg.len` bases as the UMI. Returns `(new_id, seq_off,
/// qual_off)`: the caller composes the new id (id + ":" + umi) and uses the
/// returned offset to slice seq and qual. If the read is shorter than the
/// requested UMI length, returns `None` and the read should be dropped (the
/// upstream contract is "shorter than UMI len = malformed").
#[must_use]
pub fn extract(id: &[u8], seq: &[u8], cfg: UmiConfig) -> Option<(Vec<u8>, usize)> {
    if seq.len() < cfg.len {
        return None;
    }
    let umi = &seq[..cfg.len];
    let mut new_id = Vec::with_capacity(id.len() + 1 + cfg.len);
    new_id.extend_from_slice(id);
    new_id.push(b':');
    new_id.extend_from_slice(umi);
    Some((new_id, cfg.len))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extracts_and_appends_to_id() {
        let cfg = UmiConfig {
            loc: UmiLoc::Read1,
            len: 4,
        };
        let (new_id, off) = extract(b"read_001", b"ACGTNNNN", cfg).expect("ok");
        assert_eq!(new_id, b"read_001:ACGT");
        assert_eq!(off, 4);
    }

    #[test]
    fn rejects_read_shorter_than_umi() {
        let cfg = UmiConfig {
            loc: UmiLoc::Read1,
            len: 8,
        };
        assert!(extract(b"x", b"ACGT", cfg).is_none());
    }
}
