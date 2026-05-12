# Single-cell analysis core

> Normalization, dimensionality reduction, neighbor graph construction,
> clustering, marker detection, and doublet detection — the analytical
> layer between a raw count matrix and downstream biology.

## Scope

The Scanpy / Seurat / SingleCellExperiment workflow. Covers:

- Normalization and HVG selection.
- PCA, t-SNE, UMAP.
- k-NN graph + Leiden/Louvain clustering.
- Marker gene detection / cluster annotation.
- Doublet detection (Scrublet, DoubletFinder, scDblFinder).

Spatial-specific analysis is in [`spatial.md`](spatial.md). Trajectory
methods are in [`trajectory.md`](trajectory.md). Batch correction lives
in [`integration.md`](integration.md). Multimodal joint analysis is in
[`multiomics.md`](multiomics.md).

## Design notes

- This is the layer where a strong **Rust scientific-computing stack
  pays off**: `ndarray` + `nalgebra` for linear algebra, `linfa` for
  classical ML (k-means, GMM, PCA), `petgraph` for the k-NN graph,
  `annembed` for UMAP-class embeddings, and `polars` for cells/feature
  metadata.
- The largest open piece is a **Scanpy-equivalent façade** — an
  `rsomics-sc` crate that exposes the standard pipeline
  (`pp.normalize_total`, `pp.highly_variable_genes`, `tl.pca`,
  `pp.neighbors`, `tl.leiden`, `tl.umap`, `tl.rank_genes_groups`) over
  an `anndata-rs` AnnData. SnapATAC2 already does an equivalent thing
  for scATAC; we can borrow patterns directly.
- `SingleRust` (early-stage) and `SnapATAC2` (mature, scATAC-focused)
  show that the Python+Rust pyo3 split is the realistic delivery model —
  Rust internals, Python user surface — until a Rust-native notebook
  story matures.
- Doublet detection is mostly per-cell scoring on simulated multiplets,
  which parallelises trivially with `rayon`. A pure-Rust Scrublet
  replacement is a small, isolated win.
- Leiden has multiple Rust implementations of varying maturity (none yet
  is the obvious canonical choice). Audit them before picking.

## TODO

- [~] **`Scanpy`** — Python toolkit defining the conventional scRNA
  analysis workflow (normalize → HVG → PCA → neighbors → Leiden → UMAP →
  markers).
  - Reference impl: `Python` · [scverse/scanpy](https://github.com/scverse/scanpy) · `BSD-3-Clause`
  - Existing Rust: partial — [`SingleRust`](https://github.com/SingleRust/SingleRust)
    is an early-stage Rust toolkit aiming at scanpy parity;
    [`SnapATAC2`](https://github.com/scverse/SnapATAC2) shows the Rust+pyo3
    pattern for the same domain (scATAC).
  - Existing non-C alternatives: none.
  - Priority: `P0`
  - Notes: The largest single rsomics deliverable in the single-cell
    space is a `rsomics-sc` crate that wraps `anndata-rs` + `linfa` +
    `petgraph` + `annembed` + a Leiden crate into a Scanpy-equivalent
    API. Either rebuild from scratch or fork `SingleRust`.

- [ ] **`Seurat`** — R toolkit, the other dominant single-cell
  framework.
  - Reference impl: `R / C++` · [satijalab/seurat](https://github.com/satijalab/seurat) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: none.
  - Priority: `P1`
  - Notes: Seurat's contributions (SCTransform, IntegrateData, FindMarkers)
    are well-documented. Practical rsomics path: produce a `rsomics-sc`
    that can read/write `.h5Seurat` and `.rds` Seurat objects through
    `extendr`, not a full rewrite.

- [ ] **`SingleCellExperiment` / `Bioconductor` ecosystem** — R class
  + dozens of analysis packages (scran, scater, etc.).
  - Reference impl: `R` · [Bioconductor SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) · `GPL-3+`
  - Existing Rust: none direct. AnnData ↔ SingleCellExperiment
    conversion is well-trodden via `anndata2ri` (Python) and the
    new `anndataR` package (R).
  - Existing non-C alternatives: none.
  - Priority: `P1`
  - Notes: Don't try to port — focus on round-trip IO between
    `anndata-rs` and SCE via `anndataR`.

- [~] **PCA / SVD for sparse single-cell matrices**
  - Reference impl: `scikit-learn` truncated SVD, `irlba` (R).
  - Existing Rust: [`linfa-reduction`](https://crates.io/crates/linfa-reduction)
    PCA via LOBPCG (handles high-dimensional data efficiently).
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Audit `linfa-reduction` on sparse `ndarray` / `nalgebra-sparse`
    matrices at 1M+ cells × 30K genes. May need a randomized SVD path
    (Halko-Martinsson-Tropp) for the largest atlases.

- [~] **k-NN graph (`pp.neighbors` equivalent)**
  - Reference impl: PyNNDescent / FAISS / HNSW (Python / C++).
  - Existing Rust: HNSW via [`hnsw_rs`](https://crates.io/crates/hnsw_rs),
    used by `annembed`. PyNNDescent has no direct Rust port.
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: HNSW is the right default. Wire into the rsomics-sc neighbors
    primitive, with `petgraph` storage for downstream graph ops.

- [ ] **`Leiden` clustering** — modularity-optimising community detection.
  - Reference impl: `C++ / Python` · [vtraag/leidenalg](https://github.com/vtraag/leidenalg) · `GPL-3.0`
  - Existing Rust: multiple early-stage crates (e.g. `leiden-rs`,
    `rustleiden`). None obviously canonical as of 2026-05.
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Audit existing Rust Leiden crates against `leidenalg` on
    benchmark graphs; fork the best and harden, rather than starting
    from scratch. `petgraph` integration is the priority.

- [~] **`Louvain` clustering** — predecessor to Leiden, still widely cited.
  - Reference impl: `C++ / Python` · python-louvain · `BSD-3-Clause`
  - Existing Rust: graph-clustering crates exist (`graph-clustering`,
    parts of `linfa-clustering`); maturity varies.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Lower priority than Leiden — Leiden subsumes most use cases.
    Worth supporting for legacy reproducibility only.

- [~] **`UMAP`** — non-linear embedding.
  - Reference impl: `Python` · [lmcinnes/umap](https://github.com/lmcinnes/umap) · `BSD-3-Clause`
  - Existing Rust: [`annembed`](https://crates.io/crates/annembed) — pure-Rust,
    HNSW-based embedder that produces UMAP-comparable output and is
    10× faster on large data per its paper.
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Adopt `annembed` as the default UMAP-equivalent. A `umap-rs`
    crate exists but is less mature.

- [ ] **`t-SNE`** — older non-linear embedding.
  - Reference impl: `C++` · [lvdmaaten/bhtsne](https://github.com/lvdmaaten/bhtsne) · `BSD-4-Clause`
  - Existing Rust: a small crate `bhtsne` (rust port) exists.
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: UMAP / annembed are preferred for scRNA. Keep t-SNE for
    legacy figures only.

- [ ] **Marker detection (`rank_genes_groups`)**
  - Reference impl: Scanpy / Seurat (Wilcoxon, t-test, logreg).
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Tiny pure-Rust win. `statrs` + `polars` covers the
    statistics; `rayon` parallelises across genes. Match Scanpy's
    p-value and rank ordering exactly.

- [ ] **`Scrublet`** — simulated-doublet KNN-based doublet scorer.
  - Reference impl: `Python` · [AllonKleinLab/scrublet](https://github.com/AllonKleinLab/scrublet) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Small algorithm, naturally parallel; obvious Rust rewrite
    target once the rsomics-sc pipeline exists.

- [ ] **`DoubletFinder`** — Seurat-side doublet detector.
  - Reference impl: `R` · [chris-mcginnis-ucsf/DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) · `CC BY-NC 4.0`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Non-commercial license blocks a derivative crate. Wrap via
    `extendr` for users staying in the Seurat ecosystem.

- [ ] **`scDblFinder`** — modern Bioconductor doublet detector.
  - Reference impl: `R` · [Bioconductor scDblFinder](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Better-performing than DoubletFinder in benchmarks. Use via
    `extendr`; rewrite is low priority.
