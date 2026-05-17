# rsomics-fastqc

Per-file FASTQ quality-control report — an independent Rust
reimplementation of [FastQC]. Emits `fastqc_data.txt` + `summary.txt`
(the [MultiQC]-parseable layout) per input.

## Install

```
cargo install rsomics-fastqc
```

Single binary. Auto-handles `.fq`, `.fastq`, `.fq.gz`, `.fq.bz2`,
`.fq.xz`, `.fq.zst` via [needletail].

## Usage

```
rsomics-fastqc -o qc/ reads.fastq.gz      # writes qc/reads.fastq.gz_fastqc/
rsomics-fastqc --stdout reads.fq | ...    # fastqc_data.txt to stdout
```

## Why one crate (not one per module)

The per-function partition rule asks "what does a user invoke on its
own". A FastQC *report* is that unit — nobody runs "per-base N content"
in isolation; they run the QC report and read its modules together, and
MultiQC consumes the single `fastqc_data.txt`. So, like "view a BAM",
this is one coherent operation that functionally spans modules. Per the
repo's crate-partition rule the module list and this bundling rationale
are documented explicitly here.

### Modules (FastQC contract)

1. Basic Statistics — counts, length range, %GC, encoding (never warn/fail)
2. Per base sequence quality — warn LQ<10 / median<25; fail LQ<5 / median<20
3. Per sequence quality scores — warn mode<27; fail mode<20
4. Per base sequence content — warn |A−T|/|G−C|>10%; fail >20%
5. Per sequence GC content — warn Σ|dev| >15% reads; fail >30%
6. Per base N content — warn any pos >5%; fail >20%
7. Sequence length distribution — warn unequal lengths; fail any zero-length
8. Sequence duplication levels — warn non-unique >20%; fail >50%
9. Overrepresented sequences — warn any >0.1%; fail any >1%
10. Adapter content — warn any adapter >5%; fail >10%
11. Kmer content — 7-mer binomial; warn p<0.01; fail p<1e-5
12. Per tile sequence quality — warn tile >2 below; fail >5 below

All 12 modules are implemented with the exact FastQC pass/warn/fail
thresholds. A few module *display layouts* (the long-read base-group
binning of per-position modules; the GC theoretical curve; the exact
binomial p-values; the full 209-entry contaminant label table) are
calibrated against the FastQC binary at the black-box compat gate before
the crate is published — the pass/warn/fail decision for every module is
the documented contract and is exact independent of that calibration.

## Origin

This crate is an independent Rust reimplementation of `FastQC` based on:

- The public FastQC documentation (the per-module analysis definitions
  and the exact pass/warn/fail thresholds) and the public
  `fastqc_data.txt` / `summary.txt` output format.
- Black-box behaviour comparison against the upstream FastQC binary.

FastQC is GPL-3.0 (Java); **no FastQC source code was read or used** —
the thresholds and formats come from its public Help documentation and
observed black-box output. Architectural ideas (the per-module dispatch
shape) were informed by [RastQC] (MIT), which is cited rather than
forked. Test fixtures are independently generated.

License: MIT OR Apache-2.0. Upstream credit: [FastQC] (GPL-3.0),
[RastQC] (MIT).

### External-dep quadrant classification

- `needletail` — Quadrant ① (pure Rust + SIMD).
- `rsomics-common`, `rsomics-help`, `clap`, `serde`, `serde_json`,
  `anyhow` — Quadrant ④ (edge utilities).

No FFI wrappers (no Quadrant ②); no known single-threaded-in-hot-path
deps (no Quadrant ③).

[FastQC]: https://github.com/s-andrews/FastQC
[RastQC]: https://github.com/Huang-lab/RastQC
[MultiQC]: https://multiqc.info
[needletail]: https://crates.io/crates/needletail
