# Spliced alignment (bulk RNA-seq)

> Mapping RNA-seq reads to a genome while detecting and respecting splice
> junctions.

## Scope

Covers RNA-only / splice-aware aligners. Unspliced DNA aligners (BWA,
bowtie2) live in [`02-genomics/alignment-short-read.md`](../02-genomics/alignment-short-read.md).
Pseudo / selective alignment (Salmon, kallisto) is in
[`quantification.md`](quantification.md) because its output is abundance,
not a BAM. Long-read RNA aligners (minimap2 `-x splice`) are listed but the
core minimap2 entry is in module 02.

## Design notes

- STAR's two-pass algorithm and on-disk genome suffix array are the
  performance baseline; any rewrite needs SIMD-accelerated seed extension
  and careful memory layout — the 30 GB index footprint is the single
  most painful aspect of STAR in pipelines.
- HISAT2's hierarchical FM-index is more memory-frugal and a more
  Rust-friendly target (Rust has solid FM-index crates in `rust-bio`).
  A pure-Rust HISAT-style aligner is the most tractable splice-aware
  rewrite.
- TopHat2 / MapSplice are legacy and only ported for archival / pipeline
  compatibility reasons. Authors of TopHat (Trapnell) actively redirect
  users to HISAT2.
- Long-read splice mapping is already well served by `minimap2-rs` (FFI).
  A pure-Rust long-read splice aligner would be a separate engineering
  effort and is lower priority than a short-read solution.

## TODO

- [ ] **`STAR`** — splice-aware short-read aligner with built-in chimeric
  read detection.
  - Reference impl: `C++` · [alexdobin/STAR](https://github.com/alexdobin/STAR) · `MIT`
  - Existing Rust: none pure-Rust. STAR-Fusion bindings exist only as
    pipeline wrappers.
  - Existing non-C alternatives: STARsolo (extension, same codebase, C++).
  - Priority: `P0`
  - Notes: Memory-hungry but the de-facto pipeline standard
    (ENCODE, GTEx, TCGA). Pure-Rust rewrite is a Phase-2+ project; an FFI
    binding is the realistic first step. Output is a CellRanger-compatible
    BAM + SJ.out.tab — both already covered by `noodles`.

- [ ] **`HISAT2`** — hierarchical FM-index spliced aligner.
  - Reference impl: `C++` · [DaehwanKimLab/hisat2](https://github.com/DaehwanKimLab/hisat2) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: HISAT-genotype (same authors, same code).
  - Priority: `P0`
  - Notes: GPL constrains derivative crates — a clean-room Rust port can
    be MIT/Apache. FM-index primitives in `rust-bio` and `awry`-style
    sampled suffix arrays are a credible foundation. Probably the most
    realistic full short-read splice aligner to attempt.

- [ ] **`TopHat2`** — first widely-used RNA-seq spliced aligner; built on
  bowtie2.
  - Reference impl: `C++` · [DaehwanKimLab/tophat](https://github.com/DaehwanKimLab/tophat) · `Artistic-2.0 / BSL-style`
  - Existing Rust: none.
  - Existing non-C alternatives: HISAT2 is the official successor.
  - Priority: `P2`
  - Notes: Officially deprecated by its authors (Kim, Pertea, Salzberg)
    in favour of HISAT2. Listed only because legacy GEO / SRA submissions
    still reference TopHat2 BAMs; no port needed.

- [ ] **`Subjunc`** — junction-aware aligner from the Subread package.
  - Reference impl: `C` · [Subread / Rsubread](https://subread.sourceforge.net/) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: none.
  - Priority: `P2`
  - Notes: Used mostly by R users via Rsubread; standalone usage is rare.
    Worth at least an `extendr` wrapper rather than a rewrite — the C
    code is small and stable.

- [ ] **`MapSplice`** — junction-aware aligner used heavily by TCGA legacy
  pipelines.
  - Reference impl: `C++ / Perl` · [MapSplice2 mirror](https://github.com/davidroberson/MapSplice2) · `unspecified open source`
  - Existing Rust: none.
  - Existing non-C alternatives: none.
  - Priority: `P2`
  - Notes: Largely superseded by STAR + HISAT2. Keep for TCGA replication
    work; not a porting target.
