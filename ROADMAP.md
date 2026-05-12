# Roadmap

A sketch of *order* across modules. Within each module, the per-topic TODOs
already carry P0/P1/P2 — this file is about which modules depend on which.

## Dependency picture

```
                ┌──────────────────────────────┐
                │  01-foundations              │
                │  (IO, compression, indexing, │
                │   k-mer hashing, parallelism)│
                └──────────────┬───────────────┘
                               │ used by everything below
        ┌──────────────────────┼──────────────────────┐
        │                      │                      │
        ▼                      ▼                      ▼
 02-genomics            03-transcriptomics       06-metagenomics
 (DNA alignment,        (RNA alignment,          (classification,
  variant calling,       quantification,          MAG assembly)
  assembly)              DE, splicing)
        │                      │                      │
        │                      ▼                      │
        │              04-single-cell                 │
        │              (scRNA, scATAC,                │
        │               spatial, multiomics)          │
        │                                             │
        ▼                                             │
 05-epigenomics                                       │
 (peak calling,                                       │
  methylation, Hi-C)                                  │
                                                      │
        ┌──────────────────────┬──────────────────────┘
        ▼                      ▼
 07-proteomics-         08-phylogenetics-
 structure              popgen

                  09-workflow-utility (cross-cutting)
```

## Phases

**Phase 0 — Foundations (must come first).** Solid IO and indexing crates.
Most of this exists in `noodles` already; the work is auditing and filling
gaps. Without this, every downstream tool re-invents BAM parsing.

**Phase 1 — Genomics short-read core.** Short-read aligner + variant caller +
preprocessing (fastp equivalent). Highest impact for the broadest user base.
Lots of mature C to displace.

**Phase 2 — Transcriptomics + single-cell counting.** RNA aligners and
quantifiers, then `alevin-fry`-style single-cell counting. `alevin-fry` is
already Rust — we adopt rather than rewrite.

**Phase 3 — Downstream analysis (DE, clustering, trajectory).** Statistical
analysis layer. This is where Python (Scanpy, scikit-learn) dominates;
displacing it requires good interop (PyO3) rather than feature parity from
day one.

**Phase 4 — Specialty domains.** Epigenomics, metagenomics, proteomics,
phylogenetics. Each is a self-contained subgraph that can be picked up by
contributors with domain expertise.

**Phase 5 — Workflows and viz.** Once the libraries exist, build the
Snakemake/Nextflow-style orchestrator and viz tooling.

## Non-goals

- We are **not** building a single monolithic library. Each domain is its
  own crate, composable but independently usable.
- We are **not** rewriting tools whose Rust port is already production-grade
  (`noodles`, `alevin-fry`, `needletail`, `minimap2-rs`, `nf-core` parts).
  Those are listed in the TODOs as `[x] adopt`.
- We are **not** targeting GPU as a first-class platform initially. CPU
  performance with `rayon` + SIMD comes first; GPU offload is a Phase 4+
  concern.
