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
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: maybe — spatial-neighbor enrichments parallelise
  - Upstream license: `BSD-3-Clause`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-spatial`)
  - Consumes primitives: `anndata-rs`, `zarrs`, `petgraph`, `ndarray-stats`, `linfa`, `hnsw_rs`, future `rsomics-stats`
  - Notes: Squidpy defines the API conventions (`obsp["spatial_connectivities"]`, `tl.spatial_neighbors`, `gr.nhood_enrichment`). Build the rsomics equivalent over `anndata-rs` and match the field naming exactly so Python users can read rsomics outputs in Squidpy.

- [ ] **`Giotto`** — R toolbox for spatial transcriptomics.
  - Reference impl: `R / C++` · [drieslab/Giotto](https://github.com/drieslab/Giotto) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: Squidpy (Python)
  - Parallelism: R BiocParallel + C++ inner loops
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — same as Squidpy
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-spatial` (interop via extendr)
  - Consumes primitives: `extendr`-bridge, `anndata-rs`
  - Notes: Slower than Squidpy in benchmarks (~10×). Worth supporting via `extendr` rather than rewriting — the user base is in R.

- [ ] **`stLearn`** — spatial scRNA toolkit emphasising spatial-morphology-aware clustering.
  - Reference impl: `Python` · [BiomedicalMachineLearning/stLearn](https://github.com/BiomedicalMachineLearning/stLearn) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: maybe — image feature CNNs are GPU-friendly
  - Upstream license: `BSD-3-Clause`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-spatial`
  - Consumes primitives: `anndata-rs`, image-feature backbone via `candle`/`burn`, `linfa`
  - Notes: Image-feature-augmented analysis is a niche use case. Lower priority than Squidpy parity.

- [~] **`SpaceRanger`** — 10x Visium pipeline.
  - Reference impl: `Rust + Python` · 10x Genomics · restricted
  - Existing Rust: internals are Rust at 10x but not redistributable
  - Existing Rust kind: `partial-port` (via the open `scan-rs` subset)
  - Existing non-C alternatives: open community pipelines wrapping STAR + spot-level UMI dedup
  - Parallelism: rayon + Python orchestration
  - SIMD: auto-vectorize
  - Quadrant: —
  - GPU-amenable: maybe — same constraints as Cell Ranger
  - Upstream license: restricted
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-spaceranger` as output-compatible wrapper around `alevin-fry`)
  - Consumes primitives: `alevin-fry`, `simpleaf`, `anndata-rs`
  - Notes: Match outputs (`filtered_feature_bc_matrix.h5`, `spatial/tissue_positions.csv`) from an open rsomics pipeline similar in spirit to `alevin-fry` for scRNA.

- [ ] **`Stereopy`** — STOmics' spatial transcriptomics toolkit, with particular support for Stereo-seq.
  - Reference impl: `Python` · [STOmics/Stereopy](https://github.com/STOmics/Stereopy) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: Squidpy (less Stereo-seq specific)
  - Parallelism: Python multiprocessing
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: maybe — multi-sample comparison parallelises
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-spatial`
  - Consumes primitives: `anndata-rs`, `zarrs`, `petgraph`, `ndarray-stats`
  - Notes: Important for Stereo-seq adoption. The multi-sample comparison workflows are the distinguishing feature.

- [ ] **`SpatialData`** — scverse data model for multimodal spatial experiments.
  - Reference impl: `Python` · [scverse/spatialdata](https://github.com/scverse/spatialdata) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — data model layer
  - Upstream license: `BSD-3-Clause`
  - Priority: `P0`
  - Layer: `A` (foundation — `rsomics-spatialdata` mirroring the SpatialData on-disk layout)
  - Consumes primitives: `zarrs`, `anndata-rs`, image crates
  - Notes: Adopt the SpatialData on-disk Zarr layout in `rsomics-spatial`. `zarrs` (Rust Zarr crate) covers the IO; the typed model (`Image`, `Labels`, `Points`, `Shapes`) maps cleanly to Rust enums.

- [~] **`Visium HD` tooling** — high-resolution Visium analysis (2 µm bins).
  - Reference impl: SpaceRanger 3.x + Scanpy / Seurat extensions
  - Existing Rust: none specific
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: depends on backend
  - SIMD: depends on backend
  - Quadrant: —
  - GPU-amenable: maybe — large matrices benefit from GPU sparse linear algebra
  - Upstream license: depends on backend
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-spatial`
  - Consumes primitives: `anndata-rs` out-of-core HDF5, `zarrs`
  - Notes: Mostly a question of binning resolution and very large matrices. `anndata-rs`'s out-of-core HDF5 backing matters here.

- [ ] **`Baysor`** — Bayesian segmentation for imaging-based spatial transcriptomics (MERFISH, seqFISH, Xenium).
  - Reference impl: `Julia` · [kharchenkolab/Baysor](https://github.com/kharchenkolab/Baysor) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: Baysor itself (Julia) is the only non-C entry; BOMS, RNA2seg are Python
  - Parallelism: Julia threading
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — Bayesian per-molecule assignment is SIMT-friendly per molecule
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-baysor` or fold as a segmenter inside `rsomics-spatial`)
  - Consumes primitives: `nalgebra` (multivariate Gaussians), `rayon`, `anndata-rs`, future `rsomics-stats`
  - Notes: **Best Rust target in the spatial sub-area.** The Bayesian per-molecule cell assignment naturally parallelises with `rayon` and the Gaussian / multinomial likelihoods are small `nalgebra` workloads. A Rust Baysor replacement would slot directly into SpatialData pipelines.

- [ ] **MERFISH / seqFISH / Xenium analysis tooling** — vendor pipelines (Vizgen post-processing, 10x Xenium) plus open community tools.
  - Reference impl: mixed (Python / Java / C++)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: vendor-specific
  - SIMD: vendor-specific
  - Quadrant: —
  - GPU-amenable: maybe — image-segmentation upstream is GPU-friendly
  - Upstream license: mixed (vendor-specific)
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-spatial` (vendor-specific input readers)
  - Consumes primitives: future `rsomics-baysor`, `anndata-rs`, `zarrs`
  - Notes: The unified entry point is the `rsomics-spatial` crate with a Baysor-equivalent segmenter; downstream is the same as Squidpy over an `AnnData`.
