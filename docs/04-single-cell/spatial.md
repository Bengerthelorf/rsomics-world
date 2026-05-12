# Spatial transcriptomics

> Analysis of transcripts measured with spatial coordinates: sequencing-
> based (Visium, Stereo-seq, Slide-seq) and imaging-based (MERFISH,
> seqFISH, Xenium).

## Scope

Upstream (FASTQ / image → spot/cell matrix) and downstream (spatial-
aware clustering, niche detection, segmentation, deconvolution). Bulk
spatial deconvolution into cell types is included; pure-RNA scRNA
analysis is in [`analysis-core.md`](analysis-core.md).

## Design notes

- Spatial is the youngest big sub-area of single-cell and the most
  promising place for a Rust entrant. Existing tools are split: Squidpy
  / Stereopy / stLearn / SpatialData are Python; Giotto is R; Baysor
  is Julia.
- The natural rsomics deliverable is `rsomics-spatial`, an
  `anndata-rs`-backed crate using the `SpatialData` data model (zarr +
  AnnData + image), with spatial-neighbor primitives in `petgraph` and
  spatial statistics (Moran's I, Geary's C, neighborhood enrichment) in
  `ndarray-stats`.
- Cell segmentation on imaging-based assays (Baysor, BOMS, RNA2seg) is
  algorithmically rich (Bayesian per-molecule assignment) and a great
  Rust target — `rayon` parallelises over molecules and `nalgebra`
  handles the multivariate Gaussians.
- 10x SpaceRanger is closed-source-ish; we treat it the same way as
  Cell Ranger — emit matching outputs from an open Rust pipeline.

## TODO

- [ ] **`Squidpy`** — scalable spatial omics analysis framework.
  - Reference impl: `Python` · [scverse/squidpy](https://github.com/scverse/squidpy) · `BSD-3-Clause`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Squidpy defines the API conventions (`obsp["spatial_connectivities"]`,
    `tl.spatial_neighbors`, `gr.nhood_enrichment`). Build the rsomics
    equivalent over `anndata-rs` and match the field naming exactly so
    Python users can read rsomics outputs in Squidpy.

- [ ] **`Giotto`** — R toolbox for spatial transcriptomics.
  - Reference impl: `R / C++` · [drieslab/Giotto](https://github.com/drieslab/Giotto) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: Squidpy (Python).
  - Priority: `P1`
  - Notes: Slower than Squidpy in benchmarks (~10×). Worth supporting
    via `extendr` rather than rewriting — the user base is in R.

- [ ] **`stLearn`** — spatial scRNA toolkit emphasising
  spatial-morphology-aware clustering.
  - Reference impl: `Python` · [BiomedicalMachineLearning/stLearn](https://github.com/BiomedicalMachineLearning/stLearn) · `BSD-3-Clause`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Image-feature-augmented analysis is a niche use case. Lower
    priority than Squidpy parity.

- [~] **`SpaceRanger`** — 10x Visium pipeline.
  - Reference impl: `Rust + Python` · 10x Genomics · `restricted`
  - Existing Rust: internals are Rust at 10x but not redistributable.
  - Existing non-C alternatives: open community pipelines wrapping STAR
    + spot-level UMI dedup.
  - Priority: `P1`
  - Notes: Match outputs (`filtered_feature_bc_matrix.h5`,
    `spatial/tissue_positions.csv`) from an open rsomics pipeline
    similar in spirit to `alevin-fry` for scRNA.

- [ ] **`Stereopy`** — STOmics' spatial transcriptomics toolkit, with
  particular support for Stereo-seq.
  - Reference impl: `Python` · [STOmics/Stereopy](https://github.com/STOmics/Stereopy) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: Squidpy (less Stereo-seq specific).
  - Priority: `P1`
  - Notes: Important for Stereo-seq adoption. The multi-sample comparison
    workflows are the distinguishing feature.

- [ ] **`SpatialData`** — scverse data model for multimodal spatial
  experiments.
  - Reference impl: `Python` · [scverse/spatialdata](https://github.com/scverse/spatialdata) · `BSD-3-Clause`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Adopt the SpatialData on-disk Zarr layout in `rsomics-spatial`.
    `zarrs` (Rust Zarr crate) covers the IO; the typed model
    (`Image`, `Labels`, `Points`, `Shapes`) maps cleanly to Rust enums.

- [~] **`Visium HD` tooling** — high-resolution Visium analysis (2 µm bins).
  - Reference impl: SpaceRanger 3.x + Scanpy / Seurat extensions.
  - Existing Rust: none specific.
  - Priority: `P1`
  - Notes: Mostly a question of binning resolution and very large
    matrices. `anndata-rs`'s out-of-core HDF5 backing matters here.

- [ ] **`Baysor`** — Bayesian segmentation for imaging-based spatial
  transcriptomics (MERFISH, seqFISH, Xenium).
  - Reference impl: `Julia` · [kharchenkolab/Baysor](https://github.com/kharchenkolab/Baysor) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: Baysor itself (Julia) is the only
    non-C entry; BOMS, RNA2seg are Python.
  - Priority: `P0`
  - Notes: **Best Rust target in the spatial sub-area.** The Bayesian
    per-molecule cell assignment naturally parallelises with `rayon`
    and the Gaussian / multinomial likelihoods are small `nalgebra`
    workloads. A Rust Baysor replacement would slot directly into
    SpatialData pipelines.

- [ ] **MERFISH / seqFISH / Xenium analysis tooling** — vendor pipelines
  (Vizgen post-processing, 10x Xenium) plus open community tools.
  - Reference impl: mixed (Python / Java / C++).
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: The unified entry point is the `rsomics-spatial` crate with
    a Baysor-equivalent segmenter; downstream is the same as Squidpy
    over an `AnnData`.
