# Single-cell multi-omics integration

> Joint analysis of multiple modalities measured on the same cell or
> matched cells (RNA + ATAC, RNA + protein/CITE-seq, RNA + methylation),
> plus reference-mapping.

## Scope

- Paired multi-omic integration (RNA + ATAC from 10x Multiome, CITE-seq
  RNA + ADT).
- Mosaic / partial integration where modalities are present in different
  cells (MultiVI, totalVI).
- Factor models (MOFA, MOFA+).
- Within-Seurat weighted nearest-neighbour multimodal integration (WNN).
- scATAC-side companion analysis (Signac).

Pure batch correction within a single modality is in
[`integration.md`](integration.md). Spatial-omics multimodal analysis is
in [`spatial.md`](spatial.md).

## Design notes

- This sub-area is *entirely* Python (scvi-tools) and R (Seurat, Signac,
  MOFA2). No Rust presence.
- The factor models (MOFA, MOFA+) are amenable to Rust rewrites because
  the variational inference involves only standard linear algebra (no
  neural nets). `ndarray-linalg` + `statrs` covers it.
- Deep-generative models (totalVI, MultiVI) need `candle` / `burn`;
  Phase-4 work. Near term: pyo3 bridges using `anndata-rs` /
  `MuData`-rs.
- scATAC processing already has a strong Rust foothold via SnapATAC2;
  Signac is the R counterpart and is more tightly coupled with Seurat.
- MuData (the multi-modal sibling of AnnData) needs a Rust crate. There
  is an experimental `mudata-rs` in flux; this is a foundational gap.

## TODO

- [ ] **`MOFA`** — Multi-Omics Factor Analysis (original, batch).
  - Reference impl: `Python / R` · [bioFAM/MOFA](https://github.com/bioFAM/MOFA) · `LGPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: R BiocParallel
  - SIMD: BLAS
  - Quadrant: —
  - GPU-amenable: no — superseded
  - Upstream license: `LGPL-3.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Superseded by MOFA+; keep for reproducibility of older studies. No porting target.

- [ ] **`MOFA+`** — variational Bayesian factor model for multi-modal single-cell.
  - Reference impl: `Python / R` · [bioFAM/MOFA2](https://github.com/bioFAM/MOFA2) · `LGPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: R BiocParallel + GPU optional
  - SIMD: BLAS
  - Quadrant: —
  - GPU-amenable: yes — VI updates are dense linear algebra blocks
  - Upstream license: `LGPL-3.0`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-mofa`)
  - Consumes primitives: `ndarray-linalg`, `statrs`, future `rsomics-mudata` (MuData rust crate), `anndata-rs`
  - Notes: Clean variational inference, no deep-learning component. Realistic pure-Rust target — VI updates are ~standard linear-algebra blocks. Output factors slot into `obsm["X_mofa"]` of an AnnData / MuData.

- [ ] **`WNN` (Seurat v4+)** — weighted nearest-neighbour multimodal integration.
  - Reference impl: `R / C++` · [satijalab/seurat](https://github.com/satijalab/seurat) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: R BiocParallel + C++ inner loops
  - SIMD: BLAS
  - Quadrant: —
  - GPU-amenable: maybe — k-NN search is SIMT-friendly
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-sc` (multimodal neighbour primitive inside the rsomics-sc umbrella)
  - Consumes primitives: `hnsw_rs`, `petgraph`, future `rsomics-mudata`
  - Notes: Algorithmically small once a k-NN primitive exists — WNN is a per-cell weighted combination of modality-specific neighbour graphs. Pairs naturally with the `rsomics-sc` neighbours layer.

- [ ] **`totalVI`** — VAE for joint RNA + ADT (CITE-seq).
  - Reference impl: `Python` · [scverse/scvi-tools](https://github.com/scverse/scvi-tools) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: PyTorch GPU
  - SIMD: PyTorch kernels
  - Quadrant: —
  - GPU-amenable: yes — VAE training is dense DL
  - Upstream license: `BSD-3-Clause`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-integrate` (DL family, see [`integration.md`](integration.md))
  - Consumes primitives: `candle` or `burn`, future `rsomics-mudata`, future `rsomics-stats` (NB / NB-mixture decoders)
  - Notes: Deep-learning model. PyO3 bridge first; pure Rust is Phase-4 once `candle` covers the negative-binomial-plus-NB-mixture decoders.

- [ ] **`MultiVI`** — VAE for joint / mosaic RNA + ATAC.
  - Reference impl: `Python` · [scverse/scvi-tools](https://github.com/scverse/scvi-tools) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: PyTorch GPU
  - SIMD: PyTorch kernels
  - Quadrant: —
  - GPU-amenable: yes — VAE training is dense DL
  - Upstream license: `BSD-3-Clause`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-integrate`
  - Consumes primitives: same as totalVI
  - Notes: Same constraints as totalVI. PyO3 bridge near-term.

- [ ] **`Symphony`** — reference projection (cross-listed with `integration.md`).
  - Reference impl: `R / Python` · [immunogenomics/symphony](https://github.com/immunogenomics/symphony) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: R BiocParallel
  - SIMD: BLAS
  - Quadrant: —
  - GPU-amenable: maybe — see [`integration.md`](integration.md)
  - Upstream license: `GPL-3.0`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-integrate` (canonical entry in [`integration.md`](integration.md); this row is a cross-reference)
  - Consumes primitives: see [`integration.md`](integration.md) canonical entry
  - Notes: See `integration.md`. Builds on top of Harmony; pair them in one `rsomics-integrate` crate. **Cross-reference only — canonical entry is in `integration.md`.**

- [ ] **`Signac`** — Seurat extension for scATAC and multiome analysis.
  - Reference impl: `R / C++` · [stuart-lab/signac](https://github.com/stuart-lab/signac) · `MIT`
  - Existing Rust: none directly. [`SnapATAC2`](https://github.com/scverse/SnapATAC2) is the Python+Rust functional equivalent on the Python side
  - Existing Rust kind: `none` (Signac itself); SnapATAC2 covers the Rust side
  - Existing non-C alternatives: SnapATAC2 (Rust + Python)
  - Parallelism: R BiocParallel + C++ inner loops
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: no — fragment IO + interval ops, latency-bound
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-sc` (interop via extendr; Rust-side runs via SnapATAC2)
  - Consumes primitives: `extendr`-bridge, `anndata-rs`
  - Notes: For users staying in Seurat / R, ship an `extendr` bridge that surfaces SnapATAC2 / `rsomics-sc` outputs as Signac-compatible fragment + assay objects. No need to rewrite Signac itself.

- [ ] **MuData layer** — multimodal AnnData-like container.
  - Reference impl: `Python` · [scverse/mudata](https://github.com/scverse/mudata) · `BSD-3-Clause`
  - Existing Rust: experimental — `anndata-rs` ecosystem has early MuData support; specification is still moving
  - Existing Rust kind: `partial-port`
  - Existing non-C alternatives: —
  - Parallelism: inherits from the underlying HDF5 IO
  - SIMD: inherits from HDF5 codec layer (blosc SIMD)
  - Quadrant: ②
  - GPU-amenable: no — same as AnnData IO
  - Upstream license: `BSD-3-Clause`
  - Priority: `P0`
  - Layer: `A` (foundation — `rsomics-mudata` once stable; same FFI quadrant story as `rsomics-anndata`)
  - Consumes primitives: `anndata-rs`, HDF5 FFI deps
  - Notes: Underpins every entry in this file. Make sure `anndata-rs` grows a stable MuData layer before building MOFA+ / WNN / integration-related Rust crates.
