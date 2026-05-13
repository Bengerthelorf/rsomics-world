# rsomics-world

A Cargo workspace of single-binary CLI tools that displace the C/Python/R-era
bioinformatics toolchain with **modern Rust**: fearless parallelism, explicit
SIMD, a sane installer (`cargo install rsomics-<name>`), no build-system
ceremony.

Most upstream tools are single-threaded, memory-inefficient, and written in
2005-era C or pure R. Modern multicore + SIMD + GPU resources sit idle. The
goal of `rsomics-*` is to put them to work, tool by tool, while staying
binary-compatible with upstream so existing pipelines can swap in piecewise.

## Architecture

Two layers, one workspace, one git repo. See [`CONVENTIONS.md`](CONVENTIONS.md)
for the rules; the short version:

- `crates/foundation/` — **Layer A**, library-only primitives (IO, intervals,
  k-mers, FM-index, alignment cores, stats). A crate is in A iff ≥ 2 tools
  depend on it.
- `crates/tools/` — **Layer B**, each crate is one installable binary
  (`rsomics-fastp`, `rsomics-bam`, `rsomics-bwa`, …).
- Dependency direction is **B → A → external**, enforced. A never depends on
  B; B never depends on B; sharing happens through A.

Per-domain planning lives under [`docs/`](docs/):

```
docs/
├── 00-overview/             Vision, principles, ecosystem survey, benchmarking
├── 01-foundations/          IO formats, compression, indexing, data structures
├── 02-genomics/             DNA alignment, assembly, variant calling, annotation
├── 03-transcriptomics/      Bulk RNA-seq alignment, quantification, DE, splicing
├── 04-single-cell/          scRNA, scATAC, trajectory, integration, spatial
├── 05-epigenomics/          Peak calling, methylation, Hi-C, footprinting
├── 06-metagenomics/         Classification, profiling, MAG assembly, amplicon
├── 07-proteomics-structure/ MS, structure prediction, docking
├── 08-phylogenetics-popgen/ MSA, trees, population genetics
└── 09-workflow-utility/     Workflow engines, containers, visualisation
```

Each module's `README.md` is the scope file; the topic files inside carry
TODO checklists using the entry schema in [`CONVENTIONS.md`](CONVENTIONS.md).
[`TODO.md`](TODO.md) is the flat aggregated view across modules.

## Status

Monorepo workspace, **scaffolding in place**. Foundation crates land under
`crates/foundation/`, tool binaries under `crates/tools/`. The first tool
will be installable via `cargo install rsomics-<name>` once published.

## Why Rust?

See [`docs/00-overview/motivation.md`](docs/00-overview/motivation.md). Short
version: memory safety, fearless parallelism (`rayon`), modern packaging
(`cargo`), explicit SIMD (`std::simd`), and a maturing scientific stack
(`ndarray`, `polars`, `arrow`, `candle`) make Rust a credible host language
for the next generation of bioinformatics tools — much of which is still
written in 2005-era C with hand-rolled allocators and brittle Autotools
builds.

## How to read this repo

- Start with [`docs/00-overview/`](docs/00-overview/) for context and
  principles.
- Browse module READMEs for scope.
- [`TODO.md`](TODO.md) is the cross-module checklist.
- [`CONVENTIONS.md`](CONVENTIONS.md) covers architecture, the TODO schema,
  the external-dependency quadrants, license + clean-room rules, and the
  four first-class platform targets.
