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
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Superseded by MOFA+; keep for reproducibility of older
    studies. No porting target.

- [ ] **`MOFA+`** — variational Bayesian factor model for multi-modal
  single-cell.
  - Reference impl: `Python / R` · [bioFAM/MOFA2](https://github.com/bioFAM/MOFA2) · `LGPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Clean variational inference, no deep-learning component.
    Realistic pure-Rust target — VI updates are ~standard linear-algebra
    blocks. Output factors slot into `obsm["X_mofa"]` of an AnnData /
    MuData.

- [ ] **`WNN` (Seurat v4+)** — weighted nearest-neighbour multimodal
  integration.
  - Reference impl: `R / C++` · [satijalab/seurat](https://github.com/satijalab/seurat) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Algorithmically small once a k-NN primitive exists — WNN
    is a per-cell weighted combination of modality-specific neighbour
    graphs. Pairs naturally with the `rsomics-sc` neighbours layer.

- [ ] **`totalVI`** — VAE for joint RNA + ADT (CITE-seq).
  - Reference impl: `Python` · [scverse/scvi-tools](https://github.com/scverse/scvi-tools) · `BSD-3-Clause`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Deep-learning model. PyO3 bridge first; pure Rust is Phase-4
    once `candle` covers the negative-binomial-plus-NB-mixture decoders.

- [ ] **`MultiVI`** — VAE for joint / mosaic RNA + ATAC.
  - Reference impl: `Python` · [scverse/scvi-tools](https://github.com/scverse/scvi-tools) · `BSD-3-Clause`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Same constraints as totalVI. PyO3 bridge near-term.

- [ ] **`Symphony`** — reference projection (cross-listed with
  `integration.md`).
  - Reference impl: `R / Python` · [immunogenomics/symphony](https://github.com/immunogenomics/symphony) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: See `integration.md`. Builds on top of Harmony; pair them
    in one `rsomics-integrate` crate.

- [ ] **`Signac`** — Seurat extension for scATAC and multiome analysis.
  - Reference impl: `R / C++` · [stuart-lab/signac](https://github.com/stuart-lab/signac) · `MIT`
  - Existing Rust: none directly. [`SnapATAC2`](https://github.com/scverse/SnapATAC2)
    is the Python+Rust functional equivalent on the Python side.
  - Existing non-C alternatives: SnapATAC2 (Rust + Python).
  - Priority: `P1`
  - Notes: For users staying in Seurat / R, ship an `extendr` bridge
    that surfaces SnapATAC2 / `rsomics-sc` outputs as Signac-compatible
    fragment + assay objects. No need to rewrite Signac itself.

- [ ] **MuData layer** — multimodal AnnData-like container.
  - Reference impl: `Python` · [scverse/mudata](https://github.com/scverse/mudata) · `BSD-3-Clause`
  - Existing Rust: experimental — `anndata-rs` ecosystem has early MuData
    support; specification is still moving.
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Underpins every entry in this file. Make sure `anndata-rs`
    grows a stable MuData layer before building MOFA+ / WNN /
    integration-related Rust crates.
