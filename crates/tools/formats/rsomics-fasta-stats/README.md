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

## JSON output schema (`--json`)

```jsonc
{
  "schema_version": "1.0",
  "tool": "rsomics-fasta-stats",
  "tool_version": "0.3.0",
  "status": "ok",
  "result": [                      // one element per input file
    {
      "file": "chr22.fa",          // path as supplied on CLI
      "format": "FASTA",
      "type": "DNA",               // DNA | RNA | Protein | Other
      "num_seqs": 1,               // record count
      "sum_len": 50818468,         // total bases, gaps stripped
      "min_len": 50818468,
      "max_len": 50818468,
      "avg_len": 50818468.0,       // f64, %.1f when printed
      "extended": {                // present iff --all
        "Q1": 50818468.0,          // length quartiles (float)
        "Q2": 50818468.0,
        "Q3": 50818468.0,
        "sum_gap": 0,              // bases matching --gap-letters
        "N50": 50818468,           // verbatim seqkit semantics
        "L50": 1,                  // unique-length bucket count (seqkit
                                   // tabular header still reads `N50_num`
                                   // for byte-equality compat)
        "GC(%)": 36.22,            // f64, %.2f when printed
        "sum_n": 11658691          // N (nucleotide) or X (protein) count
      }
    }
  ]
}
```

Failure envelope routes to stderr (stdout stays parseable):

```jsonc
{
  "schema_version": "1.0",
  "tool": "rsomics-fasta-stats",
  "tool_version": "0.3.0",
  "status": "error",
  "error": { "kind": "InvalidInput", "message": "..." },
  "exit_code": 1
}
```

`schema_version` is `MAJOR.MINOR`. MINOR adds optional fields, MAJOR
removes/renames. Pin against MAJOR.

## Performance

Benchmark results live in `.autopilot/state/bench-rsomics-fasta-stats-*.toml`
and the `benches/` directory. The contract for this crate: every release
must show a strictly faster wall-clock vs `seqkit stats` on the
benchmark fixtures, measured with `hyperfine --warmup 3`.
