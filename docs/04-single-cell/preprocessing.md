# Single-cell preprocessing

> From raw scRNA / scATAC / multiome FASTQs to a barcode-corrected,
> UMI-deduplicated count matrix (cells × features).

## Scope

Demultiplexing, barcode correction, UMI deduplication, transcriptome /
peak quantification at the single-cell level. Downstream analysis
(normalization, HVG, PCA, clustering) lives in
[`analysis-core.md`](analysis-core.md). Spatial-specific upstream
processing (Visium / Stereo-seq / MERFISH) is in
[`spatial.md`](spatial.md).

## Design notes

- This is the most Rust-mature corner of single-cell. `alevin-fry` /
  `simpleaf` is production-grade for scRNA, `piscem` is the recommended
  mapper, `oarfish` covers long-read scRNA. Adopt; do not rewrite.
- The two open gaps: (1) a clean Rust replacement for the gigantic
  `cellranger`/`cellranger-atac`/`cellranger-arc` Python+Rust binaries,
  and (2) a Rust scATAC barcode → fragment file pipeline that matches
  Cell Ranger ATAC output. SnapATAC2 covers downstream analysis but
  starts from a fragments file rather than FASTQ.
- 10x Genomics' own `scan-rs` (Rust) is the closest public reference for
  Cell Ranger internals — useful for behaviour matching even though it
  is not a drop-in replacement.
- For scRNA the right "pipeline default" is **piscem → alevin-fry**;
  the rsomics layer should wrap `simpleaf` rather than reimplement.

## TODO

- [~] **`Cell Ranger`** — 10x Genomics' end-to-end scRNA pipeline
  (alignment, barcode handling, gene quantification, filtering).
  - Reference impl: `Rust + Python` · [10XGenomics/cellranger](https://github.com/10XGenomics/cellranger) · `restricted (research use), some components MIT/BSD via scan-rs`
  - Existing Rust: partial — Cell Ranger itself is largely Rust
    internally; supporting library [`scan-rs`](https://github.com/10XGenomics/scan-rs)
    is open. Full pipeline is not redistributable as a Rust crate.
  - Existing non-C alternatives: `STARsolo`, `alevin-fry`,
    `kallisto|bustools`.
  - Priority: `P0`
  - Notes: Output BAM + filtered matrix h5 is a *de facto* standard.
    The realistic rsomics role is to produce indistinguishable outputs
    via `alevin-fry` / `simpleaf` so downstream Scanpy / Seurat scripts
    work unchanged.

- [ ] **`STARsolo`** — single-cell extension of STAR, drop-in Cell Ranger
  replacement.
  - Reference impl: `C++` · [alexdobin/STAR](https://github.com/alexdobin/STAR) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: alevin-fry pipeline.
  - Priority: `P1`
  - Notes: Shares all engineering pain with STAR itself (see
    `../03-transcriptomics/alignment-spliced.md`). Output is widely
    adopted in publications; matching its filtered matrix h5 exactly
    matters for benchmark reproducibility.

- [x] **`alevin-fry`** — barcode + UMI + EM-based scRNA quantification
  on top of salmon / piscem.
  - Reference impl: `Rust` · [COMBINE-lab/alevin-fry](https://github.com/COMBINE-lab/alevin-fry) · `BSD-3-Clause`
  - Existing Rust: alevin-fry (this row).
  - Existing non-C alternatives: STARsolo, Cell Ranger.
  - Priority: `P0`
  - Notes: Adopt as-is. Together with `simpleaf` and `piscem` this is
    the de-facto Rust scRNA preprocessing stack. rsomics packages it
    and contributes upstream.

- [x] **`simpleaf`** — Rust orchestrator for the alevin-fry pipeline.
  - Reference impl: `Rust` · [COMBINE-lab/simpleaf](https://github.com/COMBINE-lab/simpleaf) · `BSD-3-Clause`
  - Existing Rust: simpleaf (this row).
  - Existing non-C alternatives: none.
  - Priority: `P0`
  - Notes: Adopt. Right place to add 10x-multiome and feature-barcode
    convenience wrappers if rsomics needs them.

- [~] **`kallisto | bustools`** — pseudoalignment + barcode handling
  scRNA pipeline.
  - Reference impl: `C++` · [pachterlab/kallisto](https://github.com/pachterlab/kallisto), [BUStools/bustools](https://github.com/BUStools/bustools) · `BSD-2-Clause`
  - Existing Rust: partial — [`rust-pseudoaligner`](https://github.com/10XGenomics/rust-pseudoaligner)
    implements the kallisto T-DBG primitive used inside Cell Ranger.
  - Existing non-C alternatives: alevin-fry pipeline.
  - Priority: `P1`
  - Notes: A Rust BUS format reader/writer (`bus-rs` or similar) would
    let `rust-pseudoaligner` slot into the same downstream tooling. The
    BUS format is small and well-specified.

- [~] **`cellranger-atac`** — scATAC FASTQ → fragments + peak matrix.
  - Reference impl: `Rust + Python` · 10x Genomics · `restricted`
  - Existing Rust: most internals are Rust at 10x but not redistributable.
  - Existing non-C alternatives: `chromap` (C++, very fast),
    `SnapATAC2` (Rust, but starts from fragments).
  - Priority: `P0`
  - Notes: There is a clear opening for a Rust scATAC pipeline that goes
    from FASTQ to a 10x-compatible fragments.tsv.gz + peak/cell matrix.
    `chromap` is the fastest aligner; a Rust binary wrapping it (or a
    pure-Rust port — see `alevin-fry-atac` below) + barcode handling
    + Tn5 shift would close this gap.

- [~] **`alevin-fry-atac`** — emerging COMBINE-lab scATAC pipeline.
  - Reference impl: `Rust` · [COMBINE-lab/alevin-fry-atac](https://github.com/COMBINE-lab/alevin-fry-atac) · `BSD-3-Clause`
  - Existing Rust: alevin-fry-atac (this row, early-stage).
  - Existing non-C alternatives: cellranger-atac, chromap-based pipelines.
  - Priority: `P1`
  - Notes: Track upstream; contribute as it matures. Could become the
    canonical Rust scATAC entry point once stable.

- [~] **`cellranger-arc`** — 10x multiome (RNA + ATAC) joint pipeline.
  - Reference impl: `Rust + Python` · 10x Genomics · `restricted`
  - Existing Rust: 10x internals (not redistributable).
  - Existing non-C alternatives: none feature-complete.
  - Priority: `P1`
  - Notes: A multiome orchestrator wrapping `alevin-fry` (RNA) +
    `alevin-fry-atac` (ATAC) + a barcode-matching step is the obvious
    Rust deliverable. Output: AnnData / MuData via `anndata-rs`.
