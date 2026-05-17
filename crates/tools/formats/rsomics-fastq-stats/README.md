# rsomics-fastq-stats

Per-file statistics for FASTQ inputs. Drop-in replacement for `seqkit
stats` on FASTQ, including the quality columns (`Q20(%)`, `Q30(%)`,
`AvgQual`).

## Install

```
cargo install rsomics-fastq-stats
```

Single binary. Auto-handles `.fq`, `.fastq`, `.fq.gz`, `.fq.bz2`,
`.fq.xz`, `.fq.zst` via [needletail].

## Usage

```
rsomics-fastq-stats reads.fq
rsomics-fastq-stats --tabular --all reads.fastq.gz
rsomics-fastq-stats --json reads.fq | jq .result
```

Default columns:

```
file  format  type  num_seqs  sum_len  min_len  avg_len  max_len
```

With `--all` (verbatim seqkit FASTQ column order):

```
Q1  Q2  Q3  sum_gap  N50  N50_num  Q20(%)  Q30(%)  AvgQual  GC(%)  sum_n
```

`AvgQual` = `-10·log₁₀(Σ 10^(-q/10) / sum_len)` over every quality base
(seqkit semantics). Quality encoding defaults to sanger (offset 33);
`--fq-encoding` selects solexa / illumina-1.3+/1.5+ (offset 64). A
quality byte below the offset is a fail-loud error, not a silent clamp.

## Origin

This crate is an independent Rust reimplementation of `seqkit stats`
(FASTQ) based on:

- The seqkit paper: Shen, W. et al. *SeqKit: a cross-platform and
  ultrafast toolkit for FASTA/Q file manipulation.* PLoS ONE 11.10
  (2016) [doi:10.1371/journal.pone.0163962].
- The public FASTQ format specification.
- Black-box behaviour comparison via `--tabular` output against the
  upstream `seqkit stats` binary (the N50/L50/quartile math and the
  Q20/Q30/AvgQual formulas were locked against `seqkit v2.13.0`).

The compat contract is defined for well-formed FASTQ. On out-of-contract
degenerate input — an empty file, or zero-length reads — this crate
fails loud (no-record input) or reports zeros rather than reproducing
seqkit's `NaN` / `Unlimit` output; that is a deliberate, family-wide
choice (`rsomics-fasta-stats` behaves the same) under the fail-loud
rule: a `NaN` GC% is a wrong statistic, not a useful one.

seqkit is MIT-licensed, so clean-room is not strictly required for
licence purposes; we still document the methodology so the contract is
explicit and reproducible.

Test fixtures are independently generated; the hand-crafted tiny FASTQ
under `tests/golden/` was authored for this crate, not extracted from
seqkit's test corpus.

License: MIT OR Apache-2.0. Upstream credit: [seqkit] (MIT).

### External-dep quadrant classification

- `needletail` — Quadrant ① (pure Rust + SIMD).
- `rsomics-seqstats` — Quadrant ① (Layer-A: pure Rust length/composition
  primitives, shared with `rsomics-fasta-stats`).
- `rsomics-common`, `rsomics-help`, `clap`, `serde`, `serde_json`,
  `anyhow` — Quadrant ④ (edge utilities).

No FFI wrappers (no Quadrant ②); no known single-threaded-in-hot-path
deps (no Quadrant ③). The seqkit upstream is Go, single-threaded for
this subcommand — the perfgate provenance under `.autopilot/state/`
records the measured ratio at one thread for fairness.

[needletail]: https://crates.io/crates/needletail
[seqkit]: https://github.com/shenwei356/seqkit

## JSON output schema (`--json`)

```jsonc
{
  "schema_version": "1.0",
  "tool": "rsomics-fastq-stats",
  "tool_version": "0.1.0",
  "status": "ok",
  "result": [                      // one element per input file
    {
      "file": "reads.fq",          // path (or basename with -b)
      "format": "FASTQ",
      "type": "DNA",               // DNA | RNA | Protein | Unlimit (seqkit's
                                   // names; guessed from the first record)
      "num_seqs": 100000,
      "sum_len": 15000000,
      "min_len": 150,
      "max_len": 150,
      "avg_len": 150.0,            // f64, %.1f when printed
      "extended": {                // present iff --all
        "Q1": 150.0,               // length quartiles (float)
        "Q2": 150.0,
        "Q3": 150.0,
        "sum_gap": 0,
        "N50": 150,                // verbatim seqkit semantics
        "L50": 1,                  // unique-length bucket count (seqkit
                                   // tabular header still reads `N50_num`)
        "Q20(%)": 97.3,            // bases phred≥20 / sum_len ·100
        "Q30(%)": 93.1,            // bases phred≥30 / sum_len ·100
        "AvgQual": 35.42,          // -10·log10(Σerr / sum_len), %.2f
        "GC(%)": 50.01,            // f64, %.2f when printed
        "sum_n": 0
      }
    }
  ]
}
```

Failure envelope routes to stderr (stdout stays parseable):

```jsonc
{
  "schema_version": "1.0",
  "tool": "rsomics-fastq-stats",
  "tool_version": "0.1.0",
  "status": "error",
  "error": { "kind": "InvalidInput", "message": "..." },
  "exit_code": 1
}
```

`schema_version` is `MAJOR.MINOR`. MINOR adds optional fields, MAJOR
removes/renames. Pin against MAJOR.

## Performance

The release contract: every tagged version must show strictly faster
wall-clock than `seqkit stats -a` on the perfgate fixture (single
thread, same input). Provenance lives in `.autopilot/state/perf-*.md`
and `benches/`.
