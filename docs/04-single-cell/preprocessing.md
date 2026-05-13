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

- [~] **`Cell Ranger`** — 10x Genomics' end-to-end scRNA pipeline (alignment, barcode handling, gene quantification, filtering).
  - Reference impl: `Rust + Python` · [10XGenomics/cellranger](https://github.com/10XGenomics/cellranger) · restricted (research use); some components MIT/BSD via [`scan-rs`](https://github.com/10XGenomics/scan-rs) (binary tool, install from source — not on crates.io)
  - Existing Rust: partial — Cell Ranger itself is largely Rust internally; supporting library `scan-rs` is open. Full pipeline is not redistributable as a Rust crate.
  - Existing Rust kind: `partial-port` (the open `scan-rs` covers ~10% of the surface)
  - Existing non-C alternatives: `STARsolo`, `alevin-fry`, `kallisto|bustools`
  - Parallelism: upstream rayon + Python orchestration
  - SIMD: auto-vectorize
  - Quadrant: —
  - GPU-amenable: maybe — alignment kernel SIMT-friendly; barcode handling latency-bound
  - Upstream license: restricted (research use)
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-cellranger` as an output-compatible wrapper around `simpleaf` + `alevin-fry`)
  - Consumes primitives: `simpleaf`, `alevin-fry`, `piscem`, future `rsomics-anndata` for the filtered-matrix h5 output
  - Notes: Output BAM + filtered matrix h5 is a *de facto* standard. The realistic rsomics role is to produce indistinguishable outputs via `alevin-fry` / `simpleaf` so downstream Scanpy / Seurat scripts work unchanged. `scan-rs` is not on crates.io (squatted by `entropyscan-rs`); install from 10x source. Quadrant is `—` because the restricted upstream is not adoptable as a Rust crate; the future `rsomics-cellranger` adopts `alevin-fry`'s Quadrant ① instead.

- [ ] **`STARsolo`** — single-cell extension of STAR, drop-in Cell Ranger replacement.
  - Reference impl: `C++` · [alexdobin/STAR](https://github.com/alexdobin/STAR) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `alevin-fry` pipeline
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — alignment SIMT-friendly; cell-barcode demux less so
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-star` (single STAR binary with `--single-cell` mode rather than a separate tool)
  - Consumes primitives: same as STAR ([`../03-transcriptomics/alignment-spliced.md`](../03-transcriptomics/alignment-spliced.md))
  - Notes: Shares all engineering pain with STAR itself. Output is widely adopted in publications; matching its filtered matrix h5 exactly matters for benchmark reproducibility.

- [x] **`alevin-fry`** — barcode + UMI + EM-based scRNA quantification on top of salmon / piscem.
  - Reference impl: `Rust` · [COMBINE-lab/alevin-fry](https://github.com/COMBINE-lab/alevin-fry) · `BSD-3-Clause`
  - Existing Rust: [`alevin-fry`](https://crates.io/crates/alevin-fry) `0.14.0`
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: `STARsolo`, Cell Ranger
  - Parallelism: rayon
  - SIMD: auto-vectorize on EM loops
  - Quadrant: ①
  - GPU-amenable: maybe — EM iteration is dense; barcode resolution latency-bound
  - Upstream license: `BSD-3-Clause`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Adopt as-is. Together with `simpleaf` and `piscem` this is the de-facto Rust scRNA preprocessing stack. rsomics packages it and contributes upstream.

- [x] **`simpleaf`** — Rust orchestrator for the alevin-fry pipeline.
  - Reference impl: `Rust` · [COMBINE-lab/simpleaf](https://github.com/COMBINE-lab/simpleaf) · `BSD-3-Clause`
  - Existing Rust: [`simpleaf`](https://crates.io/crates/simpleaf) `0.24.0`
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: —
  - Parallelism: rayon + delegation to wrapped tools
  - SIMD: inherits from delegated tools
  - Quadrant: ①
  - GPU-amenable: no — orchestration layer
  - Upstream license: `BSD-3-Clause`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Adopt. Right place to add 10x-multiome and feature-barcode convenience wrappers if rsomics needs them.

- [~] **`kallisto | bustools`** — pseudoalignment + barcode handling scRNA pipeline.
  - Reference impl: `C++` · [pachterlab/kallisto](https://github.com/pachterlab/kallisto), [BUStools/bustools](https://github.com/BUStools/bustools) · `BSD-2-Clause`
  - Existing Rust: partial — [`debruijn_mapping`](https://github.com/10XGenomics/rust-pseudoaligner) (install-from-source, see [`../03-transcriptomics/quantification.md`](../03-transcriptomics/quantification.md)) implements the kallisto T-DBG primitive used inside Cell Ranger
  - Existing Rust kind: `partial-port`
  - Existing non-C alternatives: `alevin-fry` pipeline
  - Parallelism: rayon
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: maybe — T-DBG traversal irregular; barcode/UMI handling latency-bound
  - Upstream license: `BSD-2-Clause`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-kallisto-sc`, shares the kallisto core with bulk; see [`../03-transcriptomics/quantification.md`](../03-transcriptomics/quantification.md))
  - Consumes primitives: `debruijn` ecosystem, `noodles-fastq`, future `rsomics-bus` (BUS format)
  - Notes: A Rust BUS format reader/writer (`bus-rs` or similar) would let `debruijn_mapping` slot into the same downstream tooling. The BUS format is small and well-specified.

- [~] **`cellranger-atac`** — scATAC FASTQ → fragments + peak matrix.
  - Reference impl: `Rust + Python` · 10x Genomics · restricted
  - Existing Rust: most internals are Rust at 10x but not redistributable
  - Existing Rust kind: `partial-port` (via the open `scan-rs` subset)
  - Existing non-C alternatives: `chromap` (C++, very fast), `SnapATAC2` (Rust, but starts from fragments)
  - Parallelism: rayon + Python orchestration
  - SIMD: auto-vectorize; chromap's hand SIMD when wrapped
  - Quadrant: —
  - GPU-amenable: maybe — alignment SIMT-friendly; Tn5 shift CPU
  - Upstream license: restricted
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-cellranger-atac` as an output-compatible wrapper)
  - Consumes primitives: future `rsomics-chromap` or chromap FFI, `rsomics-intervals`, future `rsomics-anndata`
  - Notes: Clear opening for a Rust scATAC pipeline that goes from FASTQ to a 10x-compatible fragments.tsv.gz + peak/cell matrix. `chromap` is the fastest aligner; a Rust binary wrapping it (or a pure-Rust port) + barcode handling + Tn5 shift would close this gap.

- [~] **`alevin-fry-atac`** — emerging COMBINE-lab scATAC pipeline.
  - Reference impl: `Rust` · [COMBINE-lab/alevin-fry-atac](https://github.com/COMBINE-lab/alevin-fry-atac) · `BSD-3-Clause`
  - Existing Rust: alevin-fry-atac binary (install from source — not on crates.io as that name; `simpleaf` is the user-facing wrapper that calls it)
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: `cellranger-atac`, `chromap`-based pipelines
  - Parallelism: rayon
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: maybe — same family as `alevin-fry`
  - Upstream license: `BSD-3-Clause`
  - Priority: `P1`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Track upstream; contribute as it matures. Could become the canonical Rust scATAC entry point once stable. Drives via `simpleaf` today.

- [~] **`cellranger-arc`** — 10x multiome (RNA + ATAC) joint pipeline.
  - Reference impl: `Rust + Python` · 10x Genomics · restricted
  - Existing Rust: 10x internals (not redistributable)
  - Existing Rust kind: `partial-port` (via `scan-rs` and the COMBINE-lab tools combined)
  - Existing non-C alternatives: none feature-complete
  - Parallelism: rayon + Python orchestration
  - SIMD: auto-vectorize
  - Quadrant: —
  - GPU-amenable: maybe — same constraints as Cell Ranger
  - Upstream license: restricted
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-cellranger-arc` as a multiome orchestrator)
  - Consumes primitives: `alevin-fry` (RNA), `alevin-fry-atac` (ATAC), future `rsomics-anndata` / `rsomics-mudata`
  - Notes: A multiome orchestrator wrapping `alevin-fry` (RNA) + `alevin-fry-atac` (ATAC) + a barcode-matching step is the obvious Rust deliverable. Output: AnnData / MuData via `anndata-rs`.
