# rsomics-world

A planning and design repository for rewriting major bioinformatics software in Rust.

This repo is **docs-only** at this stage. Each module under `docs/` describes a
domain of bioinformatics, the canonical tools it contains, and the state of
existing Rust (and Zig / C / C++ / assembly) implementations. The goal is to
plan *what* to rewrite before we write a single line of Rust, so that we avoid
duplicating mature work (e.g. `noodles`, `needletail`, `minimap2-rs`) and spend
effort where Rust has the most to offer.

## Why Rust?

See [`docs/00-overview/motivation.md`](docs/00-overview/motivation.md) for the
long version. In short: memory safety, fearless parallelism, modern packaging
(`cargo`), and a maturing scientific-computing ecosystem (ndarray, polars,
arrow, candle) make Rust a credible host language for the next generation of
bioinformatics tooling — much of which is still written in 2005-era C with
hand-rolled memory management and brittle build systems.

## Repository layout

```
rsomics-world/
├── README.md              ← you are here
├── ROADMAP.md             ← cross-module sequencing and milestones
├── CONVENTIONS.md         ← naming, crate template, licensing rules
├── TODO.md                ← flat master TODO aggregated from all modules
└── docs/
    ├── 00-overview/             Vision, principles, ecosystem survey, benchmarking
    ├── 01-foundations/          IO formats, compression, indexing, data structures
    ├── 02-genomics/             DNA alignment, assembly, variant calling, annotation
    ├── 03-transcriptomics/      Bulk RNA-seq alignment, quantification, DE, splicing
    ├── 04-single-cell/          scRNA, scATAC, trajectory, integration, spatial
    ├── 05-epigenomics/          Peak calling, methylation, Hi-C, footprinting
    ├── 06-metagenomics/         Classification, profiling, MAG assembly, amplicon
    ├── 07-proteomics-structure/ MS, structure prediction, docking
    ├── 08-phylogenetics-popgen/ MSA, trees, population genetics
    └── 09-workflow-utility/     Workflow engines, containers, visualization
```

Each module directory contains a `README.md` (scope, sub-areas, design notes)
and one or more topic files with TODO checklists.

## How to read this repo

- Start with [`docs/00-overview/`](docs/00-overview/) for context and principles.
- Browse module READMEs to see what each covers.
- Use [`TODO.md`](TODO.md) for a flat checklist view across everything.
- See [`CONVENTIONS.md`](CONVENTIONS.md) before opening a PR that adds a new
  tool entry or proposes a crate name.

## Status

Planning phase. No Rust code yet — when a crate is started, it will live in a
separate repository per the conventions document, and be linked from the
relevant module doc here.
