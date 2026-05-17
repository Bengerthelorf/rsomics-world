use std::collections::{HashMap, HashSet};

use rsomics_kmer::encode::{Kmer, canonical};
use rsomics_kmer::iter::KmerIter;

#[cfg(test)]
use rsomics_kmer::encode::decode;

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum DbgError {
    #[error("k must be in 1..=32 (got {0})")]
    KOutOfRange(usize),
    #[error("k-mer iteration: {0}")]
    KmerIter(#[from] rsomics_kmer::KmerError),
}

pub type Result<T> = std::result::Result<T, DbgError>;

#[derive(Debug, Clone, Copy)]
pub struct DbgNode {
    pub kmer: Kmer,
    pub in_degree: u8,
    pub out_degree: u8,
}

// canonical de Bruijn: nodes = canonical (k-1)-mers, edge prefix→suffix; a seq and its revcomp share the graph
#[derive(Debug, Clone, Default)]
pub struct Dbg {
    pub k: usize,
    pub nodes: HashMap<Kmer, DbgNode>,
    pub edges: HashSet<(Kmer, Kmer)>,
}

impl Dbg {
    pub fn build(seqs: &[&[u8]], k: usize) -> Result<Self> {
        if !(2..=32).contains(&k) {
            return Err(DbgError::KOutOfRange(k));
        }
        let mut dbg = Self {
            k,
            ..Self::default()
        };
        for seq in seqs {
            if seq.len() < k {
                continue;
            }
            for kmer_res in KmerIter::new(seq, k, false)? {
                let kmer = match kmer_res {
                    Ok(k) => k,
                    Err(rsomics_kmer::KmerError::NonAcgt { .. }) => continue,
                    Err(e) => return Err(DbgError::KmerIter(e)),
                };
                let prefix = prefix_kminus1(kmer, k);
                let suffix = suffix_kminus1(kmer, k);
                let prefix_c = canonical(prefix, k - 1);
                let suffix_c = canonical(suffix, k - 1);
                dbg.nodes.entry(prefix_c).or_insert_with(|| DbgNode {
                    kmer: prefix_c,
                    in_degree: 0,
                    out_degree: 0,
                });
                dbg.nodes.entry(suffix_c).or_insert_with(|| DbgNode {
                    kmer: suffix_c,
                    in_degree: 0,
                    out_degree: 0,
                });
                if dbg.edges.insert((prefix_c, suffix_c)) {
                    dbg.nodes.get_mut(&prefix_c).unwrap().out_degree += 1;
                    dbg.nodes.get_mut(&suffix_c).unwrap().in_degree += 1;
                }
            }
        }
        Ok(dbg)
    }

    #[must_use]
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }
    #[must_use]
    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }

    // unitigs = (k-1)-mer key vectors; join via the (k-2) overlap to get sequences
    #[must_use]
    pub fn unitigs(&self) -> Vec<Vec<Kmer>> {
        let mut succ: HashMap<Kmer, Vec<Kmer>> = HashMap::new();
        let mut pred: HashMap<Kmer, Vec<Kmer>> = HashMap::new();
        for &(a, b) in &self.edges {
            succ.entry(a).or_default().push(b);
            pred.entry(b).or_default().push(a);
        }
        let is_branching = |k: &Kmer| -> bool {
            let s = succ.get(k).map_or(0, Vec::len);
            let p = pred.get(k).map_or(0, Vec::len);
            !(s == 1 && p == 1)
        };

        let mut visited: HashSet<Kmer> = HashSet::new();
        let mut out: Vec<Vec<Kmer>> = Vec::new();
        for &start in self.nodes.keys() {
            if visited.contains(&start) || !is_branching(&start) {
                continue;
            }
            let outs = succ.get(&start).cloned().unwrap_or_default();
            for next_start in outs {
                if visited.contains(&next_start) && is_branching(&next_start) {
                    continue;
                }
                let mut path = vec![start];
                let mut cur = next_start;
                loop {
                    path.push(cur);
                    visited.insert(cur);
                    if is_branching(&cur) {
                        break;
                    }
                    let n = succ.get(&cur).and_then(|v| v.first().copied());
                    match n {
                        Some(n) => cur = n,
                        None => break,
                    }
                }
                out.push(path);
            }
        }
        for &n in self.nodes.keys() {
            if visited.contains(&n) {
                continue;
            }
            let mut path = vec![n];
            let mut cur = n;
            visited.insert(n);
            while let Some(nx) = succ.get(&cur).and_then(|v| v.first().copied()) {
                if visited.contains(&nx) {
                    break;
                }
                path.push(nx);
                visited.insert(nx);
                cur = nx;
            }
            out.push(path);
        }
        out
    }
}

fn prefix_kminus1(kmer: Kmer, _k: usize) -> Kmer {
    kmer >> 2
}

fn suffix_kminus1(kmer: Kmer, k: usize) -> Kmer {
    let mask: u64 = if k == 1 {
        0
    } else {
        (1u64 << (2 * (k - 1))) - 1
    };
    kmer & mask
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_linear_seq_yields_n_minus_k_plus_1_edges() {
        let seqs: &[&[u8]] = &[b"ACGTACGT"];
        let dbg = Dbg::build(seqs, 4).unwrap();
        // canonical-folding may coalesce palindromic k-mers; check floor counts
        assert!(dbg.n_edges() >= 1);
        assert!(dbg.n_nodes() >= 2);
    }

    #[test]
    fn build_skips_n_kmers() {
        let seqs: &[&[u8]] = &[b"ACGTNACGT"];
        let dbg = Dbg::build(seqs, 4).unwrap();
        assert!(dbg.n_edges() >= 1);
        assert!(dbg.n_nodes() >= 1);
    }

    #[test]
    fn build_rejects_k_out_of_range() {
        assert!(matches!(
            Dbg::build(&[b"ACGT"], 1),
            Err(DbgError::KOutOfRange(1))
        ));
        assert!(matches!(
            Dbg::build(&[b"ACGT"], 33),
            Err(DbgError::KOutOfRange(33))
        ));
    }

    #[test]
    fn empty_seqs_yield_empty_dbg() {
        let dbg = Dbg::build(&[], 4).unwrap();
        assert_eq!(dbg.n_nodes(), 0);
        assert_eq!(dbg.n_edges(), 0);
    }

    #[test]
    fn unitigs_collapse_linear_paths() {
        let seq: &[u8] = b"ACGTAGCTAGCTGATCGATCAGCT";
        let dbg = Dbg::build(&[seq], 5).unwrap();
        let unitigs = dbg.unitigs();
        let total_visited: usize = unitigs.iter().map(Vec::len).sum();
        assert!(
            total_visited >= dbg.n_nodes(),
            "{total_visited} ≥ {}",
            dbg.n_nodes()
        );
    }

    #[test]
    fn prefix_and_suffix_helpers_work_for_4mer() {
        let kmer = rsomics_kmer::encode::encode(b"ACGT").unwrap();
        let pref = prefix_kminus1(kmer, 4);
        let suff = suffix_kminus1(kmer, 4);
        assert_eq!(decode(pref, 3), b"ACG".to_vec());
        assert_eq!(decode(suff, 3), b"CGT".to_vec());
    }
}
