# rsomics-fastp

Pure-Rust FASTQ preprocessor — output-compatible with [OpenGene/fastp][upstream].

Single binary, fail-fast, no defensive programming, no FFI. Reads single-end
or paired-end FASTQ (auto-detected gzip), trims adapters and poly-G tails,
filters by quality / length / N content, and emits a fastp-compatible JSON
report with per-cycle quality + content curves.

## Install

```sh
cargo install rsomics-fastp
```

## Usage

Single-end:

```sh
rsomics-fastp -i input.fastq.gz -o out.fastq.gz -j report.json
```

Paired-end:

```sh
rsomics-fastp \
  -i in_R1.fastq.gz -I in_R2.fastq.gz \
  -o out_R1.fastq.gz -O out_R2.fastq.gz \
  -j report.json
```

Disable adapter trimming, enable poly-G trimming:

```sh
rsomics-fastp -i in.fq.gz -o out.fq.gz --adapter_sequence "" --trim_poly_g
```

## Compatibility

The output FASTQ and JSON schema are designed to match `fastp` 0.23.x closely
enough that downstream tooling (MultiQC, custom dashboards) keeps working.
`tests/compat.rs` runs `rsomics-fastp` and upstream `fastp` on the same input
and asserts identical surviving-read and surviving-base counts; CI runs this
on every push.

The CLI flag set is a subset of fastp's — flags we haven't ported yet are
simply absent rather than silently ignored. Run with `--help` for the
authoritative list.

## Origin

This crate is an independent Rust reimplementation of [`fastp`][upstream]
based on:

- The published method (Chen et al., *Bioinformatics*, 2018 —
  [DOI:10.1093/bioinformatics/bty560][doi])
- The public FASTQ file format and Illumina adapter specifications
- Black-box behavior testing against the upstream binary (see
  `tests/compat.rs`)

No source code from the upstream `fastp` distribution was used as reference
during implementation. Test fixtures under `tests/golden/` are independently
synthesized; larger fixtures used for benchmarking come from public reference
datasets (HG002, 1000 Genomes).

License: MIT OR Apache-2.0.
Upstream credit: [`OpenGene/fastp`][upstream] (MIT).

[upstream]: https://github.com/OpenGene/fastp
[doi]: https://doi.org/10.1093/bioinformatics/bty560
