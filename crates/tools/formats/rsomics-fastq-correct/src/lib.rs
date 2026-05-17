// Faithful port of lh3/bfc (MIT); single-char names + signed/unsigned casts mirror source
// so the algorithm stays diff-able. Modules: kmer↔kmer.h, count↔htab.c, correct↔correct.c.
#![allow(
    clippy::many_single_char_names,
    clippy::needless_range_loop,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::cast_precision_loss
)]

mod correct;
mod count;
mod kmer;

use std::path::Path;

use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};
use rsomics_fqgz::ChunkedWriter;
use rsomics_seqio::{OwnedRecord, open_fastq};
use serde::Serialize;

use correct::correct_one;
use count::{build_table, has_unique_kmer};

const CHUNK_RECORDS: usize = 4096;

// BFC bfc_opt_init defaults; w_* weights are not CLI-exposed (BFC: "cannot be changed on command line").
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

#[derive(Debug, Default, Clone, Serialize)]
pub struct CorrectReport {
    pub reads_in: u64,
    pub reads_out: u64,
    pub reads_dropped: u64,
    pub bases_in: u64,
    pub bases_corrected: u64,
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

    pub fn run(&self, input: &Path, output: &Path) -> Result<CorrectReport> {
        if self.cfg.k < 11 || self.cfg.k > 63 || self.cfg.k.is_multiple_of(2) {
            return Err(RsomicsError::ConfigError(format!(
                "k must be odd and in 11..=63 (BFC_MAX_KMER 63), got {}",
                self.cfg.k
            )));
        }
        let ch = build_table(&[input], self.cfg)?;
        // mode computed once here (BFC bfc_correct) — per-read would be O(reads × table).
        let mode = ch.hist_mode(self.cfg.min_cov);
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
                    match correct_one(cfg, &ch, &rec, mode) {
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
    use rsomics_seqio::OwnedRecord;

    use super::{CorrectConfig, correct_one};
    use crate::correct::{EcBase, SEQ_NT4, best_island, bfc_ec_greedy_k, seq_conv};
    use crate::count::{CountTable, KmerMap};
    use crate::kmer::{BfcKmer, bfc_hash_64};

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
            map: KmerMap::default(),
        };
        let rec = OwnedRecord {
            id: b"r".to_vec(),
            seq: b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            qual: b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
        };
        assert!(correct_one(&cfg, &ch, &rec, ch.hist_mode(cfg.min_cov)).is_none());
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
            map: KmerMap::default(),
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
        let (out, _) = correct_one(&cfg, &ch, &rec, ch.hist_mode(cfg.min_cov))
            .expect("clean high-cov read corrects");
        assert!(
            out.iter().all(u8::is_ascii_uppercase),
            "a clean read must stay uppercase (no corrections): {:?}",
            String::from_utf8_lossy(&out)
        );
    }

    // BFC bfc_ec_greedy_k: best > mode/3 and 2nd < 3 → returns pos<<2|base; else -1.
    #[test]
    fn greedy_rescue_picks_the_one_supported_substitution() {
        let k = 11;
        let truth = b"ACGTCAGTTGA"; // exactly k bases
        let mut ch = CountTable {
            k,
            map: KmerMap::default(),
        };
        let mut tx = BfcKmer::NULL;
        for &b in truth {
            tx.append(k, u64::from(SEQ_NT4[b as usize]));
        }
        for _ in 0..20 {
            ch.add(&tx, true);
        }
        // 3'-most base corrupted A→T (A=0, T=3, pos k-1 = d=0); greedy must restore A.
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
            map: KmerMap::default(),
        };
        assert_eq!(bfc_ec_greedy_k(k, 0, &corrupt, &empty), -1);
    }
}
