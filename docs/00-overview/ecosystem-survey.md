# Existing Rust bioinformatics ecosystem (snapshot)

A high-level map of what already exists. Per-module docs go deeper. Verified
crates that we plan to **adopt as-is** are starred. Last reviewed
2026-05.

## Core IO and data-structure crates

- **`noodles`** ★ — pure-Rust SAM/BAM/CRAM/VCF/BCF/GFF/FASTA/FASTQ/BED. The
  reference IO library for modern Rust bioinformatics. Actively maintained
  by zaeleus.
- **`needletail`** ★ — fast FASTA/FASTQ reader with adapter detection.
- **`rust-htslib`** — thin FFI bindings to htslib. Useful when noodles
  lacks a feature; we treat it as a transitional dependency.
- **`rust-bio`** — broad utility crate (alignment algorithms, suffix arrays,
  PSSM, statistics). Older API, less actively maintained but still
  useful as a building block.
- **`bio-types`** — basic types (Interval, Strand, etc.) shared across the
  rust-bio ecosystem.
- **`gzp`**, **`niffler`** — multithreaded gzip / format-sniffing.

## Alignment and assembly

- **`minimap2-rs`** ★ — Rust bindings to minimap2 (FFI). Pure-Rust port is
  an open need.
- **`block-aligner`**, **`wfa2`** crates — SIMD/wavefront aligners.
- **`rust-spoa`** — partial-order alignment.
- No mature pure-Rust BWA / bowtie2 yet.

## Single-cell

- **`alevin-fry`** ★ — droplet-based scRNA quantification, COMBINE-lab.
  Production-grade. Adopt.
- **`oarfish`** — long-read transcript quantification.
- **`pyroe`** — Python companion; Rust internals where it counts.

## Variant calling and analysis

- **`echtvar`** — fast VCF annotation by Brent Pedersen.
- **`vcfexpress`** — VCF expression language.
- **`vartrix`** — single-cell variant tracking (10x Genomics, Rust).

## Metagenomics

- **`sourmash`** ★ — k-mer sketching for metagenomics; Rust core + Python
  bindings. Adopt.
- **`finch`** — minhash sketching.
- **`nthash`** — ntHash rolling hash.

## Workflow

- **`nf-core/rust-bioinformatics`** — community recipes.
- **`prefligth`** — pipeline preflight checks.

## Adjacent ecosystems we depend on

- **`ndarray`**, **`nalgebra`** — numerical arrays.
- **`polars`** — DataFrames; the path forward for tabular DE / DGE outputs.
- **`arrow`** — columnar interchange; aligns with Python via `pyarrow`.
- **`candle`** — pure-Rust deep learning; relevant for AlphaFold-class work.
- **`burn`** — alternative DL framework.

## What is missing (high-level)

Things below have no widely-used Rust implementation as of 2026-05:

- Production BWA / bowtie2 replacement (only FFI bindings).
- HISAT2 / STAR equivalent for spliced alignment.
- DESeq2 / edgeR — the statistics layer is firmly R.
- MACS2/3 — peak calling.
- Bismark — bisulfite alignment.
- AlphaFold inference — `candle`/`burn` make it possible, no one has
  packaged it.
- Snakemake/Nextflow equivalent — several attempts (`snakegrass`,
  `cwl-rs`) but none dominant.

Each of these gets a TODO entry under the appropriate module.
