# Trajectory inference and RNA velocity

> Pseudotime, lineage, and dynamical-modelling methods for inferring
> cell-state transitions from single-cell data.

## Scope

- Pseudotime / lineage tools that operate on a low-dimensional
  embedding (Monocle3, PAGA, Slingshot).
- RNA velocity methods that use spliced / unspliced ratios to infer
  derivatives (velocyto, scVelo).
- Aggregator frameworks (dynverse) and benchmarking suites.

Static clustering and neighbor graphs are in
[`analysis-core.md`](analysis-core.md).

## Design notes

- This is one of the most Python / R-centric corners of single-cell.
  Monocle3 (R), Slingshot (R), PAGA (Python/scanpy), scVelo (Python),
  dynverse (R). Almost nothing Rust.
- The numerical primitives needed are mostly already covered by
  `petgraph` (PAGA-style graph contraction), `linfa` (principal curves
  for Slingshot), and `ndarray-stats` (dynamic-EM for scVelo). What is
  missing is the wrapping logic and an `AnnData` ↔ Rust round-trip.
- scVelo's dynamical model uses ODE-likelihood fitting per gene — a
  good `rayon` parallelisation target if rewritten. `candle` or `burn`
  are not strictly needed.
- Pragmatic ordering: pseudotime tools wrap via `extendr` / pyo3 first;
  RNA velocity is the most promising native-Rust target because the
  preprocessing (counting spliced vs. unspliced) is heavy IO that
  `noodles-bam` does well.

## TODO

- [ ] **`Monocle3`** — graph-abstraction-based pseudotime.
  - Reference impl: `R` · [cole-trapnell-lab/monocle3](https://github.com/cole-trapnell-lab/monocle3) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Builds on UMAP + a custom graph-clustering principal-graph
    learner. The principal-graph step is the algorithmically novel
    piece; Rust port is feasible with `petgraph` + `linfa` but is
    non-trivial. First step: wrap via `extendr` and ship a Rust
    AnnData ↔ Monocle3 bridge.

- [ ] **`PAGA`** — Partition-based graph abstraction (Scanpy).
  - Reference impl: `Python` · [scverse/scanpy (paga)](https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.paga.html) · `BSD-3-Clause`
  - Existing Rust: none.
  - Notes: Algorithmically tiny — given a clustering and a neighbor
    graph, aggregate edges between clusters and test connectivity.
    `petgraph` + statistical-test crate covers it. Good early Rust
    rewrite candidate inside `rsomics-sc`.

- [ ] **`Slingshot`** — principal-curves pseudotime over a clustering.
  - Reference impl: `R` · [kstreet13/slingshot](https://github.com/kstreet13/slingshot) · `Artistic-2.0`
  - Existing Rust: none. Python port [`pyslingshot`](https://github.com/mossjacob/pyslingshot)
    exists.
  - Existing non-C alternatives: pyslingshot (Python).
  - Priority: `P1`
  - Notes: Lineage tree from cluster MST + smoothed principal curves.
    Algorithm is compact; if `linfa` gets a principal-curve primitive
    this rewrites into a few hundred lines.

- [ ] **`scVelo`** — dynamical RNA velocity model with latent time.
  - Reference impl: `Python` · [theislab/scvelo](https://github.com/theislab/scvelo) · `BSD-3-Clause`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: scVelo only supports Python by design — the EM over the
    splicing-kinetics ODE per gene is the bottleneck and an obvious
    `rayon` win. Pair with a Rust spliced/unspliced counter (see
    `velocyto` below). Long-term scVelo replacement is a strong
    Phase-3 target.

- [ ] **`velocyto`** — per-cell spliced / unspliced counting from BAM.
  - Reference impl: `Python + C++` · [velocyto-team/velocyto.py](https://github.com/velocyto-team/velocyto.py) · `BSD-2-Clause`
  - Existing Rust: none direct, but [`alevin-fry`](https://github.com/COMBINE-lab/alevin-fry)
    USA-mode quantification already emits spliced/unspliced/ambiguous
    counts per cell — the modern recommendation.
  - Existing non-C alternatives: alevin-fry USA mode (Rust).
  - Priority: `P1`
  - Notes: Mark `velocyto` itself as legacy; the rsomics pipeline
    routes through alevin-fry USA. Add a thin `polars` exporter so
    scVelo / a future Rust velocity tool can consume it directly.

- [ ] **`dynverse`** — meta-framework benchmarking 70+ trajectory
  methods.
  - Reference impl: `R` · [dynverse](https://github.com/dynverse/dynverse) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: dynverse's value is the benchmarking harness, not the
    methods. If rsomics ships its own trajectory tools, a small
    `rsomics-bench-trajectory` crate replicating dynverse's metric
    suite (TI evaluation) is a useful Phase-4 deliverable.
