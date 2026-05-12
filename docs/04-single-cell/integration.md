# Batch integration and reference mapping

> Removing technical / batch effects across samples or studies so cells
> from different batches share a comparable embedding.

## Scope

Methods that take multiple AnnData / Seurat objects (or a single AnnData
with a batch covariate) and produce a corrected embedding or a corrected
count matrix. Includes linear methods (Harmony, MNN/fastMNN), graph
methods (BBKNN), nearest-neighbour methods (Scanorama), and deep
generative models (scVI, LIGER, Symphony).

Multimodal joint integration (RNA + ATAC, CITE-seq) is in
[`multiomics.md`](multiomics.md). Spatial-aware integration is in
[`spatial.md`](spatial.md).

## Design notes

- Almost everything here is Python (Harmony / harmonypy, BBKNN,
  Scanorama, scvi-tools) or R (Seurat IntegrateData, LIGER, fastMNN).
  No mature Rust integration tooling exists.
- The linear methods (Harmony, MNN) are small enough to rewrite in pure
  Rust — Harmony in particular is essentially an iterative
  cluster-aware ridge correction over a PCA embedding and could fit in
  one focused crate using `ndarray-linalg`.
- Deep-learning methods (scVI, totalVI, MultiVI) need `candle` or `burn`.
  Reaching pure-Rust feature parity is **Phase-4** work; near-term we
  bridge via PyO3 to `scvi-tools`.
- The same `anndata-rs` AnnData representation should host both batch
  metadata and the various corrected embeddings (`obsm["X_harmony"]`,
  `obsm["X_scvi"]`, etc.), matching Scanpy convention so downstream
  tooling is unchanged.

## TODO

- [ ] **`Harmony`** — iterative PCA correction by cluster soft-assignment.
  - Reference impl: `R` · [immunogenomics/harmony](https://github.com/immunogenomics/harmony) · `GPL-3.0`
    · Python port [`slowkow/harmonypy`](https://github.com/slowkow/harmonypy)
  - Existing Rust: none.
  - Existing non-C alternatives: harmonypy (Python).
  - Priority: `P0`
  - Notes: Small algorithm (a few hundred lines of pure
    linear algebra and k-means). The single most-used scRNA batch
    correction method. Excellent first-target Rust port — entire
    algorithm fits in one `ndarray-linalg` + `linfa-clustering` crate.

- [ ] **`scVI`** — variational-autoencoder-based batch correction and
  imputation.
  - Reference impl: `Python` · [scverse/scvi-tools](https://github.com/scverse/scvi-tools) · `BSD-3-Clause`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Needs a deep-learning backend (`candle` / `burn`) and a
    `torch.distributions`-equivalent layer (negative binomial,
    zero-inflated NB). Long road. Near term: pyo3 bridge to
    scvi-tools using `anndata-rs` as the shared store.

- [ ] **`BBKNN`** — batch-balanced k-NN graph.
  - Reference impl: `Python` · [Teichlab/bbknn](https://github.com/Teichlab/bbknn) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Very small algorithm: do per-batch k-NN, splice neighbour
    lists. With `hnsw_rs` already in rsomics-sc this is a few hundred
    lines. Cheap, high-impact rewrite.

- [ ] **`Scanorama`** — panorama-stitching cell alignment.
  - Reference impl: `Python` · [brianhie/scanorama](https://github.com/brianhie/scanorama) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Algorithm is nearest-neighbour matching + SVD-based
    correction. Fits within the same `linfa` + `hnsw_rs` core that
    Harmony uses.

- [ ] **`LIGER` / `rliger`** — non-negative matrix factorisation
  integration.
  - Reference impl: `R / C++` · [welch-lab/liger](https://github.com/welch-lab/liger) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: PyLiger (Python).
  - Priority: `P2`
  - Notes: iNMF kernel is amenable to Rust (`ndarray-linalg`). Less
    commonly used than Harmony / scVI; lower priority.

- [ ] **`MNN` / `fastMNN`** — mutual nearest-neighbour batch correction.
  - Reference impl: `R / C++` · [Bioconductor batchelor](https://bioconductor.org/packages/release/bioc/html/batchelor.html) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Small algorithm; the cosine-normalised MNN search reuses the
    same HNSW infrastructure. Pair with Harmony in a single
    `rsomics-integrate` crate.

- [ ] **`Symphony`** — reference-mapping (project new cells onto an
  existing Harmony-integrated reference).
  - Reference impl: `R / Python` · [immunogenomics/symphony](https://github.com/immunogenomics/symphony) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Natural extension of a Rust Harmony — the reference-mapping
    logic is light on top of an existing Harmony embedding. Worth
    bundling once Harmony is ported.
