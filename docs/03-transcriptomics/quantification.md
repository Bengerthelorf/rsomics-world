# Transcript and gene quantification

> Turning aligned (or pseudo-aligned) reads into transcript / gene abundance
> estimates.

## Scope

Two families:

1. **Lightweight / selective-alignment quantifiers** (Salmon, kallisto,
   piscem) — index a transcriptome and produce TPM/counts directly from
   FASTQ.
2. **Counters on top of BAM** (featureCounts, HTSeq-count, RSEM) — start
   from a STAR/HISAT2 BAM and a GTF.

StringTie does both quantification and assembly; the assembly half lives
in [`assembly-isoform.md`](assembly-isoform.md).

## Design notes

- Selective alignment + EM is the area where Rust already has the most
  visible presence: COMBINE-lab's `oarfish` is pure Rust, `simpleaf` is
  pure Rust, and `piscem` is a Rust binary wrapping a C++ index. Adopting
  these and contributing back is more productive than greenfield porting.
- For BAM-based counters the bottleneck is interval-tree querying of GFF
  features. `rust-bio` has interval trees, `noodles-bam` has efficient
  iteration, so a `featureCounts` clone is a small project with high
  pipeline impact.
- HTSeq-count is Python and slow but still cited because of strict-mode
  semantics. A drop-in Rust replacement with bit-identical counts in
  intersection-strict / intersection-nonempty / union modes would be
  immediately useful.
- RSEM does EM over a transcriptome-aligned BAM. Its math is shared with
  Salmon / oarfish, but its IO layer is dated and a Rust rewrite is
  attractive.

## TODO

- [ ] **`Salmon`** — selective-alignment-based transcript quantification.
  - Reference impl: `C++` · [COMBINE-lab/salmon](https://github.com/COMBINE-lab/salmon) · `GPL-3.0`
  - Existing Rust: none directly; the successor `piscem` is a Rust binary
    around a C++ core, see below. `alevin-fry` (Rust) consumes salmon /
    piscem RAD output.
  - Existing non-C alternatives: `piscem` (Rust+C++ hybrid).
  - Priority: `P0`
  - Notes: Bulk-quant feature parity with `salmon quant` is the single
    most impactful Rust deliverable in transcriptomics. The maths
    (selective-alignment scoring + range-factor EM) is well-documented;
    the engineering pain is the compact suffix-array index. Coordinate
    with COMBINE-lab on whether `piscem` becomes the canonical path.

- [~] **`kallisto`** — pseudoalignment-based transcript quantification
  using the T-DBG.
  - Reference impl: `C++` · [pachterlab/kallisto](https://github.com/pachterlab/kallisto) · `BSD-2-Clause`
  - Existing Rust: partial — [`rust-pseudoaligner`](https://github.com/10XGenomics/rust-pseudoaligner)
    from 10x Genomics implements the T-DBG pseudoalignment primitive in
    Rust, used inside Cell Ranger but not packaged as a kallisto drop-in.
  - Existing non-C alternatives: `BUStools` (companion, C++).
  - Priority: `P0`
  - Notes: `rust-pseudoaligner` plus a EM-quantification layer would give
    most of `kallisto quant`. Output formats (BUS, h5) need `noodles` /
    custom HDF5 layers (`hdf5-rust`).

- [ ] **`RSEM`** — EM-based transcript expression estimation from
  transcriptome alignments.
  - Reference impl: `C++ / Perl` · [deweylab/RSEM](https://github.com/deweylab/RSEM) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: `salmon --alignment-mode` covers most
    use cases.
  - Priority: `P1`
  - Notes: Still mandatory for some isoform-DE pipelines that expect
    `.isoforms.results`. A Rust port would be a relatively small project
    once `noodles` + `nalgebra`/`ndarray` are in place. Match RSEM's
    log-likelihood model exactly so DESeq2 / tximport downstream is
    unaffected.

- [ ] **`featureCounts`** (Subread) — interval-based read summarization
  over a GFF.
  - Reference impl: `C` · [Subread/featureCounts](https://subread.sourceforge.net/featureCounts.html) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: `htseq-count` (Python).
  - Priority: `P0`
  - Notes: Most commonly-used gene-level counter in bulk RNA-seq.
    Small, focused tool — implement as a `noodles-bam` + interval-tree
    binary `rsomics-featurecounts`. Match output column-for-column with
    `Rsubread::featureCounts` so existing DESeq2 / edgeR scripts work
    unchanged.

- [ ] **`HTSeq-count`** — Python read-feature counter, reference for
  intersection-strict semantics.
  - Reference impl: `Python / Cython` · [htseq/htseq](https://github.com/htseq/htseq) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: `featureCounts` (faster, different
    semantics), `featurecounts-rs` would inherit.
  - Priority: `P1`
  - Notes: Slow; pipelines tolerate it for strict-mode semantics. A
    Rust drop-in matching `--mode intersection-strict` + `--mode union`
    bit-for-bit is a useful Phase-2 deliverable.

- [ ] **`StringTie`** (quantification mode) — transcript abundance over
  spliced BAMs.
  - Reference impl: `C++` · [gpertea/stringtie](https://github.com/gpertea/stringtie) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: none.
  - Priority: `P1`
  - Notes: See `assembly-isoform.md` for the assembly side; the
    `-e` (estimate only, given GTF) mode is just quantification. Same
    rewrite covers both.

- [x] **`oarfish`** — long-read transcript quantification with coverage-aware EM.
  - Reference impl: `Rust` · [COMBINE-lab/oarfish](https://github.com/COMBINE-lab/oarfish) · `BSD-3-Clause`
  - Existing Rust: oarfish (this row).
  - Existing non-C alternatives: `IsoQuant` (Python), `bambu` (R).
  - Priority: `P0`
  - Notes: Adopt as-is. COMBINE-lab uses `minimap2-rs` underneath. We
    package it in the rsomics ecosystem and contribute upstream when
    needed.

- [~] **`piscem`** — next-gen Rust-wrapped selective-alignment + index for
  Salmon / alevin-fry.
  - Reference impl: `Rust + C++` · [COMBINE-lab/piscem](https://github.com/COMBINE-lab/piscem) · `BSD-3-Clause`
  - Existing Rust: piscem (front-end Rust, internals still C++).
  - Existing non-C alternatives: salmon, kallisto.
  - Priority: `P0`
  - Notes: Adopt and track. Eventually the C++ core (compacted dBG
    index, mapping kernel) becomes the obvious port target for a fully
    pure-Rust salmon-class quantifier.
