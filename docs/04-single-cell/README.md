# 04 — Single-cell and spatial

> Single-cell RNA-seq, scATAC, multiome, and spatial transcriptomics — from
> raw FASTQ to integrated cross-modality atlases.

## Sub-topics

- [`preprocessing.md`](preprocessing.md) — Cell Ranger, STARsolo,
  alevin-fry, kallisto|bustools, cellranger-arc, cellranger-atac.
- [`analysis-core.md`](analysis-core.md) — Scanpy / Seurat /
  SingleCellExperiment: PCA, neighbor graph, Leiden/Louvain, UMAP/t-SNE,
  marker detection, doublet detection.
- [`trajectory.md`](trajectory.md) — Monocle3, PAGA, Slingshot, scVelo,
  dynverse, velocyto.
- [`integration.md`](integration.md) — Harmony, scVI, BBKNN, Scanorama,
  LIGER, MNN/fastMNN, Symphony.
- [`spatial.md`](spatial.md) — Squidpy, Giotto, stLearn, SpaceRanger,
  Stereopy, SpatialData, Baysor.
- [`multiomics.md`](multiomics.md) — MOFA / MOFA+, WNN, totalVI,
  MultiVI, Symphony, Signac.

## Cross-cutting design notes

- Single-cell is the area where Rust already has a foothold: the
  COMBINE-lab stack (`alevin-fry`, `simpleaf`, `piscem`, `oarfish`)
  for transcriptomics quantification, `SnapATAC2` for scATAC analysis
  (Rust internals + Python API), and `anndata-rs` for HDF5 backed
  on-disk AnnData. The rsomics task is to **integrate and extend**
  these, not to re-engineer from scratch.
- The analytical core (PCA, neighbor graph, Leiden, UMAP) maps to
  `linfa`, `petgraph`, `annembed`, and `linfa-clustering`. Each of those
  needs auditing and a thin scanpy-equivalent layer.
- Trajectory and integration are dominated by Python / R *deep learning*
  models (scVelo, scVI, totalVI). Rewriting them as pure-Rust models
  requires `candle` / `burn` and is **Phase-4** work; near-term, the
  interop layer (PyO3 + `anndata-rs`) is what matters.
- Spatial transcriptomics is the youngest sub-area and has the most room
  for a fresh Rust entrant — Baysor is Julia, Squidpy/Stereopy are
  Python, and Bayesian segmentation is naturally parallelizable.
- Multiome / multimodal analysis is firmly Seurat / scvi-tools territory.
  Plan for `extendr` and PyO3 bridges and a shared on-disk MuData /
  AnnData layout rather than rewrites.
