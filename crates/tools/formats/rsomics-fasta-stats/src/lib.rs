//! FASTA-only per-file statistics, byte-compatible with the FASTA subset of
//! `seqkit stats --tabular`.
//!
//! `seqkit stats` covers FASTA + FASTQ in one binary. This crate is the
//! FASTA-only partition of that surface — the FASTQ partition lives in a
//! sibling crate. The split matches the per-function partition rule that
//! governs the workspace.
//!
//! ## Reference output (seqkit, FASTA-only columns)
//!
//! Default columns: `file format type num_seqs sum_len min_len avg_len max_len`
//! `--all` adds: `Q1 Q2 Q3 sum_gap N50 N50_num GC(%) sum_n`
//!
//! FASTQ-only columns (`Q20(%)`, `Q30(%)`, `AvgQual`) are NOT emitted by this
//! crate — pass a FASTQ file to a future `rsomics-fastq-stats` instead.
//!
//! ## How compatibility is anchored
//!
//! `--tabular` output is the strict-compat surface. The default pretty form
//! mirrors seqkit's two-space-padded human output, including
//! comma-grouped integers (matching Go's `humanize.Comma`), but the compat
//! test only diffs the tabular form so that humanize edge cases don't pin
//! us to Go's exact behaviour.

pub mod compute;
pub mod output;

pub use compute::{Config, FastaStats, SeqType, compute_stats};
pub use output::{render_pretty, render_tabular};
