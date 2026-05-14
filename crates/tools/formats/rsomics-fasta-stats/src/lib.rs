//! FASTA-only per-file statistics, byte-compatible with the FASTA subset
//! of `seqkit stats --tabular`. The FASTQ partition lives in
//! `rsomics-fastq-stats`.
//!
//! Default columns: `file format type num_seqs sum_len min_len avg_len max_len`.
//! With `--all`: `Q1 Q2 Q3 sum_gap N50 N50_num GC(%) sum_n`.
//!
//! `--tabular` output is the compat anchor — the pretty form is for
//! humans, not byte-equal against seqkit.

pub mod compute;
pub mod output;

pub use compute::{Config, FastaStats, SeqType, compute_stats};
pub use output::{render_pretty, render_tabular};
