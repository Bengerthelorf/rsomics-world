# Benchmarking plan

How we measure whether a Rust rewrite is worth shipping.

## What we compare against

For every rewrite we keep two baselines:

1. **The canonical C/C++ tool** at its latest stable release (e.g. BWA-MEM
   0.7.17, samtools 1.20).
2. **The most credible existing Rust implementation**, if any (e.g. noodles
   for IO, minimap2-rs FFI for long-read alignment).

The new crate must clearly beat the Rust baseline and come within a target
% of the C baseline before we declare it production-ready. Targets vary
by category — see per-module docs.

## What we measure

- **Wall-clock time** on a fixed input on a fixed machine.
- **Peak RSS** memory.
- **CPU utilization** under thread scaling sweeps (1, 4, 16, 64, 128
  threads). Linearity matters more than single-thread perf.
- **Output equivalence** — bit-identical where possible (sorted BAM,
  normalized VCF), or content-equivalent with a published diff tool.

## Reference datasets

A small set of fixed inputs lives in a sibling `rsomics-bench-data` repo
(planned). They are public, citation-friendly, and span scales:

| Dataset | Size | Used by |
|---|---|---|
| NA12878 chr22 30x | ~1 GB | Alignment, variant calling |
| HG002 ONT chr22 50x | ~5 GB | Long-read alignment |
| 1000 Genomes phase 3 chr22 VCF | ~200 MB | VCF tools |
| PBMC 10k (10x v3) | ~25 GB | Single-cell |
| MGnify mock community | ~10 GB | Metagenomics |
| ENCODE GM12878 H3K27ac ChIP | ~2 GB | Peak calling |

Smaller smoke-test slices (~10 MB) live next to each crate for CI.

## Hardware

A canonical benchmark box is defined per-quarter. Currently:

- 64-core AMD EPYC, 256 GB RAM, NVMe scratch, no GPU.

Optional GPU rig for AlphaFold-class work, defined separately when relevant.

## CI integration

Every crate's CI runs criterion benches on smoke data. A
"nightly-benchmarks" job runs the full reference set weekly on the canonical
box and posts results to a static site. PRs that regress by >5% on the
smoke set block merge.

## Honesty rules

- We publish median, p99, and stdev — not just min.
- We publish the command lines used. Reproducibility over storytelling.
- If the C tool wins, the README says so. We do not hide losses behind
  "warm cache" runs.
