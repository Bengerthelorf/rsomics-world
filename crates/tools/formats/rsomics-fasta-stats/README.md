# rsomics-fasta-stats

Per-file statistics for FASTA inputs. Drop-in replacement for the FASTA
subset of `seqkit stats`.

## Install

```
cargo install rsomics-fasta-stats
```

Single binary. Auto-handles `.fa`, `.fasta`, `.fa.gz`, `.fa.bz2`,
`.fa.xz`, `.fa.zst` via [needletail].

## Usage

```
rsomics-fasta-stats genome.fa
rsomics-fasta-stats --tabular --all genome.fa.gz
rsomics-fasta-stats --json scaffolds.fasta | jq .result
```

Default columns:

```
file  format  type  num_seqs  sum_len  min_len  avg_len  max_len
```

With `--all`:

```
Q1  Q2  Q3  sum_gap  N50  N50_num  GC(%)  sum_n
```

FASTQ quality columns (`Q20(%)`, `Q30(%)`, `AvgQual`) are not emitted —
use `rsomics-fastq-stats` for FASTQ inputs.

## Origin

This crate is an independent Rust reimplementation of `seqkit stats`
based on:

- The seqkit paper: Shen, W. et al. *SeqKit: a cross-platform and
  ultrafast toolkit for FASTA/Q file manipulation.* PLoS ONE 11.10
  (2016) [doi:10.1371/journal.pone.0163962].
- The public FASTA format specification.
- Black-box behaviour comparison via `--tabular` output against the
  upstream `seqkit stats` binary.

seqkit is MIT-licensed, so clean-room is not strictly required for
licence purposes; we still document the methodology so the contract is
explicit and reproducible for future GPL upstreams that share the same
pipeline shape.

Test fixtures are independently generated; the hand-crafted tiny FASTA
under `tests/golden/` was authored for this crate, not extracted from
seqkit's test corpus.

License: MIT OR Apache-2.0. Upstream credit: [seqkit] (MIT).

[needletail]: https://crates.io/crates/needletail
[seqkit]: https://github.com/shenwei356/seqkit

## Performance

Benchmark results live in `.autopilot/state/bench-rsomics-fasta-stats-*.toml`
and the `benches/` directory. The contract for this crate: every release
must show a strictly faster wall-clock vs `seqkit stats` on the
benchmark fixtures, measured with `hyperfine --warmup 3`.
