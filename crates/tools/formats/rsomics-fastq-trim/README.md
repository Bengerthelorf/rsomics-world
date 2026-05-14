# rsomics-fastq-trim

Adapter, poly-G, poly-X, and fixed-length trimming for FASTQ inputs.
Single binary, gzip in/out, multi-threaded over chunks of records.

## Install

```
cargo install rsomics-fastq-trim
```

## Scope

This crate is the **trim-only** partition of fastp's surface. The
per-function partition explicitly splits fastp's bundled tasks into
sibling crates:

| Task | Crate |
|---|---|
| 3' adapter trim (static + auto-detect) | **rsomics-fastq-trim** ← here |
| PE overlap-based adapter detection | **rsomics-fastq-trim** ← here |
| Poly-G / Poly-X trim | **rsomics-fastq-trim** ← here |
| Fixed-length 5'/3' trim | **rsomics-fastq-trim** ← here |
| Quality + length + sliding-window filter | rsomics-fastq-quality |
| UMI extract / stamp | rsomics-fastq-umi |
| FASTQ summary stats | rsomics-fastq-stats |
| Adapter index k-mer detection | (future) |

Composing these is a one-liner pipeline; bundling them is the
fastp-shaped Swiss-army anti-pattern this project is unwinding.

## Usage

```
# SE adapter trim with a known fastp default
rsomics-fastq-trim -i in.fq.gz -o out.fq.gz \
    --adapter-sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

# PE overlap-detect adapter trim, 8-thread, JSON envelope
rsomics-fastq-trim -i r1.fq.gz -I r2.fq.gz -o r1.tr.fq.gz -O r2.tr.fq.gz \
    --detect-adapter-for-pe -t 8 --json | jq .result

# Poly-G trim (NovaSeq 2-color chemistry)
rsomics-fastq-trim -i in.fq.gz -o out.fq.gz --trim-poly-g

# Fixed-len trim of first 6 bp + last 4 bp
rsomics-fastq-trim -i in.fq.gz -o out.fq.gz --trim-front 6 --trim-tail 4
```

## Origin

Independent Rust reimplementation of the trim hot path from `fastp`
(MIT-licensed). Methodology:

- Source-read the upstream `adaptertrimmer.cpp`, `overlapanalysis.cpp`,
  `polyx.cpp`, and `knownadapters.h` — fastp is MIT-licensed so source
  reading is allowed and is the established practice for matching
  upstream semantics in this project.
- Compat is anchored byte-level on FASTQ output for the canonical flag
  combinations against `fastp` on the same input.
- Defaults match fastp where reasonable; deviations are documented in
  `--help` and in the JSON envelope's `tool_version`.

Test fixtures are independently generated. The chr22 paired-end fixture
used for benches lives under `tests/fixtures-staging/` (gitignored) and
is provisioned by the Tier-2 fetcher in `rsomics-common`.

License: MIT OR Apache-2.0. Upstream credit: [fastp] (MIT).

[fastp]: https://github.com/OpenGene/fastp

## JSON output schema (`--json`)

```jsonc
{
  "schema_version": "1.0",
  "tool": "rsomics-fastq-trim",
  "tool_version": "0.3.0",
  "status": "ok",
  "result": {
    "reads_in": 50000,                    // input record count (PE counts both mates)
    "reads_out": 49852,                   // emitted record count
    "bases_in": 7500000,
    "bases_out": 7421300,
    "adapter_trimmed_reads": 312,         // records where static-seq adapter fired
    "adapter_trimmed_bases": 9_240,
    "poly_g_trimmed_reads": 0,
    "poly_g_trimmed_bases": 0,
    "poly_x_trimmed_reads": 0,
    "poly_x_trimmed_bases": 0,
    "fixed_trimmed_bases": 0,             // bases removed by --trim_front* / --trim_tail*
    "overlap_trimmed_pairs": 0,           // PE pairs where overlap detect fired
    "overlap_trimmed_bases": 0,
    "reads_too_short_after_trim": 148     // post-trim length < --length_required
  }
}
```

Failure envelope routes to stderr (stdout stays parseable):

```jsonc
{
  "schema_version": "1.0",
  "tool": "rsomics-fastq-trim",
  "tool_version": "0.3.0",
  "status": "error",
  "error": { "kind": "InvalidInput", "message": "..." },
  "exit_code": 1
}
```

All counts are `u64`. `*_reads` counters count individual records (in PE
mode a single record-event on either mate counts once). `_pairs`
counters count pair-events that fired together. `schema_version` is
`MAJOR.MINOR` — pin against MAJOR.

## Performance

The per-function perf hard rule: every release must show strictly
faster wall-clock vs upstream `fastp` on the same machine, plus
multi-axis evidence (peak RSS, instructions retired, page faults).
Provenance for each release lives in
`.autopilot/state/bench-rsomics-fastq-trim-*.toml`.
