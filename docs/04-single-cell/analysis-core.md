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

- [~] **`Scanpy`** — Python toolkit defining the conventional scRNA analysis workflow (normalize → HVG → PCA → neighbors → Leiden → UMAP → markers).
  - Reference impl: `Python` · [scverse/scanpy](https://github.com/scverse/scanpy) · `BSD-3-Clause`
  - Existing Rust: partial — [`SingleRust`](https://github.com/SingleRust/SingleRust) is an early-stage Rust toolkit aiming at Scanpy parity; [`SnapATAC2`](https://github.com/scverse/SnapATAC2) shows the Rust+pyo3 pattern for the same domain (scATAC)
  - Existing Rust kind: `partial-port`
  - Existing non-C alternatives: —
  - Parallelism: rayon + ndarray on the Rust side; SnapATAC2 demonstrates pyo3 boundary
  - SIMD: auto-vectorize via ndarray BLAS
  - Quadrant: ① (in-memory ops) / ② (HDF5 IO via anndata-rs)
  - GPU-amenable: maybe — PCA / UMAP / Leiden have GPU-friendly subsets; the workflow as a whole isn't
  - Upstream license: `BSD-3-Clause`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-sc` as the Scanpy-equivalent façade)
  - Consumes primitives: `anndata-rs`, `linfa`, `linfa-reduction`, `hnsw_rs`, `petgraph`, `annembed`, `polars`, future `rsomics-stats`
  - Notes: The largest single rsomics deliverable in the single-cell space is a `rsomics-sc` crate that wraps `anndata-rs` + `linfa` + `petgraph` + `annembed` + a Leiden crate into a Scanpy-equivalent API. Either rebuild from scratch or fork `SingleRust`.

- [ ] **`Seurat`** — R toolkit, the other dominant single-cell framework.
  - Reference impl: `R / C++` · [satijalab/seurat](https://github.com/satijalab/seurat) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: R BiocParallel + C++ inner loops
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — SCTransform NB regression GPU-friendly; IntegrateData anchor finding less so
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-sc` (interop layer — read/write h5Seurat / RDS via extendr)
  - Consumes primitives: `anndata-rs`, `extendr`-bridge, `polars`
  - Notes: Seurat's contributions (SCTransform, IntegrateData, FindMarkers) are well-documented. Practical rsomics path: produce a `rsomics-sc` that can read/write `.h5Seurat` and `.rds` Seurat objects through `extendr`, not a full rewrite.

- [ ] **`SingleCellExperiment` / `Bioconductor` ecosystem** — R class + dozens of analysis packages (scran, scater, etc.).
  - Reference impl: `R` · [Bioconductor SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) · `GPL-3+`
  - Existing Rust: none direct. AnnData ↔ SingleCellExperiment conversion is well-trodden via `anndata2ri` (Python) and the new `anndataR` package (R)
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: R BiocParallel
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — the umbrella R class is metadata; specific algorithms vary
  - Upstream license: `GPL-3+`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-sc` (round-trip IO bridge)
  - Consumes primitives: `anndata-rs`, `extendr`-bridge
  - Notes: Don't try to port — focus on round-trip IO between `anndata-rs` and SCE via `anndataR`.

- [~] **PCA / SVD for sparse single-cell matrices**
  - Reference impl: `scikit-learn` truncated SVD, `irlba` (R)
  - Existing Rust: [`linfa-reduction`](https://crates.io/crates/linfa-reduction) `0.8.1` (PCA via LOBPCG, handles high-dimensional data efficiently)
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: —
  - Parallelism: rayon via ndarray-linalg
  - SIMD: BLAS-level explicit
  - Quadrant: ①
  - GPU-amenable: yes — SVD is dense linear algebra
  - Upstream license: `MIT OR Apache-2.0` (linfa)
  - Priority: `P0`
  - Layer: `adopt` (with possible randomized-SVD addition)
  - Consumes primitives: —
  - Notes: Audit `linfa-reduction` on sparse `ndarray` / `nalgebra-sparse` matrices at 1M+ cells × 30K genes. May need a randomized SVD path (Halko-Martinsson-Tropp) for the largest atlases.

- [~] **k-NN graph (`pp.neighbors` equivalent)**
  - Reference impl: PyNNDescent / FAISS / HNSW (Python / C++)
  - Existing Rust: [`hnsw_rs`](https://crates.io/crates/hnsw_rs) `0.3.4`, used by `annembed`. PyNNDescent has no direct Rust port
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: —
  - Parallelism: rayon
  - SIMD: auto-vectorize on distance computation
  - Quadrant: ①
  - GPU-amenable: yes — k-NN search is SIMT-friendly (large-scale ANN ports to GPU exist)
  - Upstream license: `MIT OR Apache-2.0`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: HNSW is the right default. Wire into the rsomics-sc neighbors primitive, with `petgraph` storage for downstream graph ops.

- [ ] **`Leiden` clustering** — modularity-optimising community detection.
  - Reference impl: `C++ / Python` · [vtraag/leidenalg](https://github.com/vtraag/leidenalg) · `GPL-3.0`
  - Existing Rust: multiple early-stage crates (e.g. `leiden-rs`, `rustleiden`). None obviously canonical as of 2026-05
  - Existing Rust kind: `partial-port`
  - Existing non-C alternatives: —
  - Parallelism: per-iteration partly parallel
  - SIMD: limited
  - Quadrant: ③ (existing pure-Rust crates are single-threaded on the hot path)
  - GPU-amenable: maybe — community-detection passes have GPU variants in the literature
  - Upstream license: `GPL-3.0`
  - Priority: `P0`
  - Layer: `A` (foundation — `rsomics-leiden` once a canonical pick emerges)
  - Consumes primitives: `petgraph`
  - Notes: Audit existing Rust Leiden crates against `leidenalg` on benchmark graphs; fork the best and harden, rather than starting from scratch. `petgraph` integration is the priority. Clean-room derivation needed because of GPL upstream.

- [~] **`Louvain` clustering** — predecessor to Leiden, still widely cited.
  - Reference impl: `C++ / Python` · python-louvain · `BSD-3-Clause`
  - Existing Rust: graph-clustering crates exist (`graph-clustering`, parts of `linfa-clustering`); maturity varies
  - Existing Rust kind: `partial-port`
  - Existing non-C alternatives: —
  - Parallelism: limited
  - SIMD: none
  - Quadrant: ③
  - GPU-amenable: maybe — same family as Leiden
  - Upstream license: `BSD-3-Clause`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-leiden` (Leiden umbrella crate with `--algorithm louvain` flag)
  - Consumes primitives: `petgraph`
  - Notes: Lower priority than Leiden — Leiden subsumes most use cases. Worth supporting for legacy reproducibility only.

- [~] **`UMAP`** — non-linear embedding.
  - Reference impl: `Python` · [lmcinnes/umap](https://github.com/lmcinnes/umap) · `BSD-3-Clause`
  - Existing Rust: [`annembed`](https://crates.io/crates/annembed) `0.1.6` — pure-Rust, HNSW-based embedder that produces UMAP-comparable output and is ~10× faster on large data per its paper; also placeholder [`umap`](https://crates.io/crates/umap) `0.1.0` ("TBD")
  - Existing Rust kind: `rust-native` (annembed's algorithm is the crate's own contribution, related-to-but-not-a-port-of UMAP)
  - Existing non-C alternatives: —
  - Parallelism: rayon
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: maybe — force-directed embedding GPU variants exist (cuML UMAP)
  - Upstream license: `MIT OR Apache-2.0`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Adopt `annembed` as the default UMAP-equivalent. The `umap` crate at 0.1.0 is a placeholder.

- [ ] **`t-SNE`** — older non-linear embedding.
  - Reference impl: `C++` · [lvdmaaten/bhtsne](https://github.com/lvdmaaten/bhtsne) · `BSD-4-Clause`
  - Existing Rust: [`bhtsne`](https://crates.io/crates/bhtsne) `0.5.4` (Rust port)
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: —
  - Parallelism: rayon
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: maybe — Barnes-Hut tree traversal irregular; quad-tree GPU variants exist
  - Upstream license: `BSD-4-Clause`
  - Priority: `P2`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: UMAP / annembed are preferred for scRNA. Keep t-SNE for legacy figures only.

- [ ] **Marker detection (`rank_genes_groups`)**
  - Reference impl: Scanpy / Seurat (Wilcoxon, t-test, logreg)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: rayon over genes
  - SIMD: auto-vectorize
  - Quadrant: —
  - GPU-amenable: no — per-gene rank stats, memory-latency-bound
  - Upstream license: `BSD-3-Clause` (Scanpy)
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-sc` (a `markers` subcommand inside the rsomics-sc umbrella)
  - Consumes primitives: `polars`, `statrs`, `rayon`, future `rsomics-stats`
  - Notes: Tiny pure-Rust win. `statrs` + `polars` covers the statistics; `rayon` parallelises across genes. Match Scanpy's p-value and rank ordering exactly.

- [ ] **`Scrublet`** — simulated-doublet KNN-based doublet scorer.
  - Reference impl: `Python` · [AllonKleinLab/scrublet](https://github.com/AllonKleinLab/scrublet) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: rayon over simulated doublets
  - SIMD: auto-vectorize
  - Quadrant: —
  - GPU-amenable: maybe — simulated-doublet generation is trivially parallel; k-NN scoring is GPU-friendly
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-sc`
  - Consumes primitives: `hnsw_rs`, `rayon`, `linfa`, future `rsomics-stats`
  - Notes: Small algorithm, naturally parallel; obvious Rust rewrite target once the rsomics-sc pipeline exists.

- [ ] **`DoubletFinder`** — Seurat-side doublet detector.
  - Reference impl: `R` · [chris-mcginnis-ucsf/DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) · `CC BY-NC 4.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: R BiocParallel
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: no — small algorithm; license blocks Rust derivative
  - Upstream license: `CC BY-NC 4.0` (non-commercial)
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Non-commercial license blocks a derivative crate. Wrap via `extendr` for users staying in the Seurat ecosystem.

- [ ] **`scDblFinder`** — modern Bioconductor doublet detector.
  - Reference impl: `R` · [Bioconductor scDblFinder](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: R BiocParallel
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: no
  - Upstream license: `GPL-3.0`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-sc` (interop via extendr)
  - Consumes primitives: `extendr`-bridge
  - Notes: Better-performing than DoubletFinder in benchmarks. Use via `extendr`; rewrite is low priority.
