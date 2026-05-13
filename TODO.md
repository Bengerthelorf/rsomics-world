# Master TODO

A flat aggregated view across modules. The authoritative checklists live in
each module's per-topic file under `docs/`. This file is a one-page index +
status dashboard.

Legend (4 values, see [`CONVENTIONS.md`](CONVENTIONS.md) for full definitions):

- `[ ]` **open** — no Rust implementation, on our queue.
- `[~]` **partial** — Rust crate exists but incomplete (FFI-wrapper, partial-port, or rust-native research-grade).
- `[x]` **adopt** — production-grade pure-Rust crate exists; adopt as a direct dependency.
- `[A]` **subprocess-adopt** — upstream is fine as-is; we invoke it as a process but don't link or rewrite.

Priorities: P0 must-have · P1 high value · P2 nice to have.

## Status dashboard (post-Phase-1 reclassification, 2026-05-14)

| Module | `[ ]` open | `[~]` partial | `[x]` adopt | `[A]` subprocess | Total |
|---|---:|---:|---:|---:|---:|
| [01 — Foundations](docs/01-foundations/) | 7 | 12 | 26 | 0 | 45 |
| [02 — Genomics](docs/02-genomics/) | 55 | 5 | 8 | 0 | 68 |
| [03 — Transcriptomics](docs/03-transcriptomics/) | 28 | 2 | 2 | 0 | 32 |
| [04 — Single-cell & spatial](docs/04-single-cell/) | 37 | 12 | 2 | 0 | 51 |
| [05 — Epigenomics](docs/05-epigenomics/) | 37 | 1 | 1 | 0 | 39 |
| [06 — Metagenomics](docs/06-metagenomics/) | 37 | 0 | 2 | 0 | 39 |
| [07 — Proteomics & structure](docs/07-proteomics-structure/) | 36 | 0 | 3 | 1 | 40 |
| [08 — Phylogenetics & popgen](docs/08-phylogenetics-popgen/) | 28 | 0 | 0 | 1 | 29 |
| [09 — Workflow & utility](docs/09-workflow-utility/) | 28 | 0 | 3 | 0 | 31 |
| **Total** | **293** | **32** | **47** | **2** | **374** |

47 tools already have production-grade pure-Rust implementations — we adopt those rather than rewrite. 32 have partial coverage (FFI wrappers like `rust-htslib`, partial ports like `ruzstd`, or rust-native research-grade tools like `rust-mdbg`). 2 we adopt via subprocess only (UShER, Foldseek). 293 are genuine open work.

### Layer assignment breakdown

| Layer | Count |
|---|---:|
| A (foundation crate) | 14 |
| B (standalone tool binary) | 115 |
| `adopt` (existing Rust crate as direct dep) | 60 |
| `subcommand-of-<crate>` (folded into an umbrella binary) | 107 |
| `—` (not pursued — legacy / GUI / proprietary / commercial / closed-source / deprecated; or cross-reference to canonical entry) | 78 |
| **Total** | **374** |

## Module entry points

### 01 — Foundations *(IO, compression, indexing, data structures, parallelism)*

- [`io-formats.md`](docs/01-foundations/io-formats.md)
- [`compression.md`](docs/01-foundations/compression.md)
- [`indexing.md`](docs/01-foundations/indexing.md)
- [`data-structures.md`](docs/01-foundations/data-structures.md)
- [`parallelism.md`](docs/01-foundations/parallelism.md)

Foundations is the most "done" module — `noodles`, `needletail`,
`flate2`, `niffler`, `rayon`, `nthash`, `sourmash`, `finch`,
`probabilistic-collections` etc. cover most of it. Real open gaps:
AnnData/h5ad full read+write, MAF (TCGA flavour) parsing, parallel
suffix-array constructor for plant-genome FM-index builds.

### 02 — Genomics *(DNA alignment, assembly, variant calling, annotation)*

- [`alignment-short-read.md`](docs/02-genomics/alignment-short-read.md)
- [`alignment-long-read.md`](docs/02-genomics/alignment-long-read.md)
- [`assembly.md`](docs/02-genomics/assembly.md)
- [`variant-calling.md`](docs/02-genomics/variant-calling.md)
- [`sv-calling.md`](docs/02-genomics/sv-calling.md)
- [`annotation.md`](docs/02-genomics/annotation.md)
- [`preprocessing.md`](docs/02-genomics/preprocessing.md)

Largest open surface. **P0 anchors:** pure-Rust BWA-MEM2 / Strobealign,
fastp equivalent, GATK HaplotypeCaller class caller, VEP-class
annotator. Use `varlociraptor` and `echtvar` as proof Rust can ship
here.

### 03 — Transcriptomics *(bulk RNA-seq)*

- [`alignment-spliced.md`](docs/03-transcriptomics/alignment-spliced.md)
- [`quantification.md`](docs/03-transcriptomics/quantification.md)
- [`assembly-isoform.md`](docs/03-transcriptomics/assembly-isoform.md)
- [`differential-expression.md`](docs/03-transcriptomics/differential-expression.md)
- [`splicing.md`](docs/03-transcriptomics/splicing.md)

The DE layer is firmly R/Python; the path forward is interop
(`extendr`, PyO3) + ndarray/polars statistics rather than pure
rewrite. **P0 anchors:** STAR / HISAT2 replacement; featureCounts
replacement.

### 04 — Single-cell & spatial

- [`preprocessing.md`](docs/04-single-cell/preprocessing.md)
- [`analysis-core.md`](docs/04-single-cell/analysis-core.md)
- [`trajectory.md`](docs/04-single-cell/trajectory.md)
- [`integration.md`](docs/04-single-cell/integration.md)
- [`spatial.md`](docs/04-single-cell/spatial.md)
- [`multiomics.md`](docs/04-single-cell/multiomics.md)

`alevin-fry`, `simpleaf`, `modkit`, `oarfish` already ship — adopt.
Real Rust opportunity: a Scanpy-equivalent `rsomics-sc` with
`annembed` UMAP, `linfa` clustering, `anndata-rs` IO. **P0 anchors:**
the Scanpy-equivalent itself; Harmony port; Baysor (already mature
Julia — port priority is open).

### 05 — Epigenomics

- [`peak-calling.md`](docs/05-epigenomics/peak-calling.md)
- [`chip-atac-pipelines.md`](docs/05-epigenomics/chip-atac-pipelines.md)
- [`methylation.md`](docs/05-epigenomics/methylation.md)
- [`chromatin-3d.md`](docs/05-epigenomics/chromatin-3d.md)
- [`footprinting.md`](docs/05-epigenomics/footprinting.md)

`modkit` is the standout existing Rust. **P0 anchors:** MACS3, cooler
(`.cool` IO + computation), MethylDackel, pairtools.

### 06 — Metagenomics

- [`classification.md`](docs/06-metagenomics/classification.md)
- [`profiling.md`](docs/06-metagenomics/profiling.md)
- [`assembly-mag.md`](docs/06-metagenomics/assembly-mag.md)
- [`amplicon.md`](docs/06-metagenomics/amplicon.md)

`sourmash`, `skani` already shipping. **P0 anchors:** Kraken2-class
classifier, MEGAHIT-class assembler, MetaBAT2-class binner.

### 07 — Proteomics & structure

- [`mass-spectrometry.md`](docs/07-proteomics-structure/mass-spectrometry.md)
- [`structure-prediction.md`](docs/07-proteomics-structure/structure-prediction.md)
- [`structure-analysis.md`](docs/07-proteomics-structure/structure-analysis.md)
- [`docking.md`](docs/07-proteomics-structure/docking.md)

`Sage` (MSFragger-class search), `mzdata` (MS IO), `pdbtbx` (PDB/mmCIF)
already exist. AlphaFold-class work needs `candle`/`burn` more than
new code. **P0 anchors:** OpenMS-class library; AlphaFold inference
pipeline; AutoDock Vina port.

### 08 — Phylogenetics & population genetics

- [`alignment-msa.md`](docs/08-phylogenetics-popgen/alignment-msa.md)
- [`trees.md`](docs/08-phylogenetics-popgen/trees.md)
- [`population-genetics.md`](docs/08-phylogenetics-popgen/population-genetics.md)

`UShER` is the existing Rust standout. **P0 anchors:** MAFFT / MUSCLE5
class MSA; IQ-TREE / RAxML class tree builder; PLINK2-class popgen
toolkit.

### 09 — Workflow & utility

- [`workflow-engines.md`](docs/09-workflow-utility/workflow-engines.md)
- [`containers.md`](docs/09-workflow-utility/containers.md)
- [`data-viz.md`](docs/09-workflow-utility/data-viz.md)

`Pixi` and `rattler` (Conda re-implementation in Rust) ship today.
`bigtools` for bigWig/bigBed too. **P0 anchors:** a Snakemake/Nextflow
equivalent in Rust; browser-side genome browser components.

---

## How to use this file

- The dashboard counts above are hand-maintained. When you add or change a
  TODO entry in a module doc, also bump the relevant row.
- Sort priority within a module by reading that module's file. There is
  no global priority ordering — different downstream pipelines need
  different P0s.
- See [`ROADMAP.md`](ROADMAP.md) for cross-module phasing.
- See [`CONVENTIONS.md`](CONVENTIONS.md) for the entry format every TODO
  must follow.
