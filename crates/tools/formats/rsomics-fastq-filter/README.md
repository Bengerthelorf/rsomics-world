# rsomics-fastq-filter

Per-read quality + length filter for FASTQ inputs. Reads pass or fail as a whole
unit — no trimming is performed. SE and PE modes (PE: both mates must pass or the
pair is dropped).

## Install

```
cargo install rsomics-fastq-filter
```

## Scope

This crate is the **filter-only** partition of fastp's surface. The
per-function partition explicitly separates fastp's bundled tasks into
sibling crates:

| Task | Crate |
|---|---|
| Per-read quality + length filter (pass/fail whole read) | **rsomics-fastq-filter** ← here |
| 3' adapter + poly-G/X + fixed-length trim | rsomics-fastq-trim |
| Quality + sliding-window cut | rsomics-fastq-quality |
| UMI extract / stamp | rsomics-fastq-umi |
| FASTQ summary stats | rsomics-fastq-stats |

## Usage

```
# SE filter with fastp defaults (Q15 threshold, 40% low-qual limit, min length 15)
rsomics-fastq-filter -i in.fq.gz -o out.fq.gz

# Strict: Q20 threshold, 20% low-qual limit, min length 50 bp
rsomics-fastq-filter -i in.fq.gz -o out.fq.gz --qualified_quality_phred 20 --unqualified_percent_limit 20 -l 50

# PE: pair is dropped if either mate fails
rsomics-fastq-filter -i r1.fq.gz -I r2.fq.gz -o r1.filt.fq.gz -O r2.filt.fq.gz

# Length-only: disable quality filter
rsomics-fastq-filter -i in.fq.gz -o out.fq.gz -Q -l 30

# JSON report envelope
rsomics-fastq-filter -i in.fq.gz -o out.fq.gz --json | jq .result
```

## Filter semantics

Sourced from `fastp` `src/filter.cpp` and `src/options.h` (fastp MIT; reading
and citing is the established practice in this project).

A base at position `i` is **unqualified** iff:
```
qual_byte[i] < (qualified_quality_phred + 33)
```
The comparison is strict less-than on the raw ASCII quality byte (fastp
`simd.cpp`, `CountQualityMetricsImpl`).

A read **fails quality** iff:
```
low_qual_count > unqualified_percent_limit * len / 100.0   (strict >)
OR
n_count > n_base_limit                                      (strict >)
```

A read **fails length** iff:
```
len < required_length                   (strict <)
OR (max_length > 0 AND len > max_length) (strict >)
```

For PE, the pair is kept only if **both** mates pass all enabled filters.

### Defaults

| Flag | Default | fastp default |
|---|---|---|
| `--qualified_quality_phred` | 15 | Q15 (`qualifiedQual='0'` ASCII) |
| `--unqualified_percent_limit` | 40 | 40 |
| `--n_base_limit` | 5 | 5 |
| `-l / --length_required` | 15 | 15 |
| `--length_limit` | 0 (no upper bound) | 0 |

## JSON output schema (`--json`)

```jsonc
{
  "schema_version": "1.0",
  "tool": "rsomics-fastq-filter",
  "tool_version": "0.1.0",
  "status": "ok",
  "result": {
    "mode": "SE",
    "reads_in": 50000,
    "reads_out": 46200,
    "bases_in": 7500000,
    "bases_out": 6930000,
    "reads_failed_quality": 2100,
    "reads_failed_length": 1400,
    "reads_failed_n_bases": 300
  }
}
```

In PE mode, each mate's failure reason is attributed individually (matching
fastp's per-read counting). When exactly one mate of a pair fails, 2 reads
are removed from `reads_out` but only 1 per-reason counter increments; the
invariant `reads_in - reads_out == sum(per-reason counters)` holds only when
both mates of every failing pair share the same failure reason.

## Origin

Independent Rust reimplementation of fastp's quality/length filter, informed by:

- The fastp paper: Chen, S. et al. *fastp: an ultra-fast all-in-one
  FASTQ preprocessor.* Bioinformatics 34.17 (2018), i884–i890
  [doi:10.1093/bioinformatics/bty560].
- Reading the upstream source: `src/filter.cpp` (`Filter::passFilter`),
  `src/options.h` (`QualityFilteringOptions`, `ReadLengthFilteringOptions`),
  and `src/simd.cpp` (`CountQualityMetricsImpl`) — fastp is MIT-licensed so
  source reading is allowed and is the established practice for matching
  upstream semantics in this project.
- Black-box behaviour comparison via byte-equal output against `fastp`
  (see `tests/compat.rs`).

License: MIT OR Apache-2.0.
Upstream credit: [fastp](https://github.com/OpenGene/fastp) (MIT).

### External-dep quadrant classification

- `rsomics-seqio` — the workspace FASTQ reader: decode-only producer thread
  + parallel parse; ISA-L igzip gz backend on Linux (Quadrant ② via the
  isolated `rsomics-igzip`), pure-Rust flate2 elsewhere (Quadrant ①).
- `flate2` (zlib-rs backend) — Quadrant ① (pure Rust + SIMD; `parallel_gz`
  output tests).
- `libdeflater` — Quadrant ② (FFI wrapper around the libdeflate C library;
  used in the chunked parallel gzip output pipeline — the hot compression
  path). Documented as the one FFI dep in this crate.
- `rayon` — Quadrant ① (pure Rust parallelism; drives the chunked parallel
  filter + compression pipeline).
- `rsomics-common`, `rsomics-help`, `clap`, `serde`, `serde_json`, `anyhow` — Quadrant ④.
