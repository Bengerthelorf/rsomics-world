# Amplicon analysis (16S / ITS / 18S)

> Per-sample ASV/OTU inference, clustering, taxonomy assignment, and
> downstream community analysis from amplicon (16S rRNA, ITS, 18S) data.

## Scope

Includes: end-to-end amplicon pipelines (QIIME2, mothur), denoisers
producing ASVs (DADA2, Deblur, UNOISE3), dereplication and OTU clustering
(VSEARCH, USEARCH, swarm), and amplicon-based functional prediction
(PICRUSt2). Excludes: shotgun community profiling (see
[classification](classification.md), [profiling](profiling.md)).

## Design notes

- The DADA2 algorithm is essentially Poisson error modeling + greedy
  partitioning of unique sequences. A pure-Rust DADA2 is a high-value,
  scoped target: tight math, no exotic dependencies, and current users
  are stuck on R for the entire downstream chain.
- VSEARCH is the open-source workhorse — it covers dereplication, OTU
  clustering, chimera detection, and pairwise alignment. The C++ code is
  ~80k lines but well-modularized. Rewrite in pieces (`rsomics-derep`,
  `rsomics-chimera`).
- QIIME2 is a Python plugin framework, not a tool — the right Rust strategy
  is to ship our amplicon tools as QIIME2-compatible plugins (or PyO3
  bindings) rather than rebuild the framework.
- Reference taxonomy assignment (SILVA, GTDB, UNITE) is a separate problem
  that overlaps with classification (Kraken2, Centrifuge). Don't duplicate.
- PICRUSt2 = phylogenetic placement (EPA-NG / SEPP) + hidden-state
  prediction (HSP). The hidden-state-prediction step is a clean Rust
  rewrite (`linfa`); the placement step inherits whatever we do for GTDB-Tk.
- License watch: QIIME2 BSD-3, DADA2 LGPL-3 (R package), mothur GPL-3,
  VSEARCH dual BSD-2/GPL-3, USEARCH proprietary, swarm GPL-3, Deblur BSD,
  PICRUSt2 GPL-3.

## TODO

- [ ] **`QIIME2`** — plugin-based amplicon analysis framework.
  - Reference impl: `Python` · [qiime2 on GitHub](https://github.com/qiime2) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing + plugin-specific
  - SIMD: inherited from plugins
  - Quadrant: —
  - GPU-amenable: no — plugin framework
  - Upstream license: `BSD-3-Clause`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: future rsomics-* binaries (as QIIME2 plugins)
  - Notes: Reproducing the framework is not the goal — interop is. Plan to ship `rsomics-dada2` (and others) as installable QIIME2 plugins, so users keep the Python provenance + visualization layer while running Rust kernels underneath.

- [ ] **`DADA2`** — sample-resolved ASV inference with Poisson error modeling.
  - Reference impl: `R` + `C++` · [benjjneb/dada2](https://github.com/benjjneb/dada2) · `LGPL-3` (R package)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `Deblur` (Python, different algorithm)
  - Parallelism: R BiocParallel + C++ inner loops
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — error-rate EM is dense; partitioning irregular
  - Upstream license: `LGPL-3`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-dada2`)
  - Consumes primitives: `noodles-fastq`, `polars`, `rayon`, future `rsomics-stats` (Poisson EM)
  - Notes: The most-used amplicon denoiser. Algorithm is well-described and Rust-friendly: error-rate EM + partitioning. Free us from the R toolchain and let users run DADA2 from a static binary or QIIME2 plugin. Largest user-facing win in this sub-area.

- [ ] **`mothur`** — classic OTU-clustering amplicon pipeline.
  - Reference impl: `C++` · [mothur/mothur](https://github.com/mothur/mothur) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: no — pre-ASV era pipeline
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Pre-ASV era pipeline. Still cited; declining usage. Skip in favor of building first-class DADA2/VSEARCH equivalents.

- [ ] **`USEARCH`** — proprietary OTU/UNOISE3 pipeline.
  - Reference impl: `C++` · drive5.com · proprietary (free 32-bit, paid 64-bit)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `VSEARCH` (open-source clone)
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: no — closed source
  - Upstream license: proprietary
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Closed-source, license-restricted. VSEARCH was created precisely to replace it. The community has largely migrated. Listed only because UNOISE3 (the denoising algorithm) is still occasionally needed; that algorithm is reproducible from the published paper.

- [ ] **`VSEARCH`** — open-source replacement for USEARCH.
  - Reference impl: `C++` · [torognes/vsearch](https://github.com/torognes/vsearch) · `BSD-2-Clause OR GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE/AVX
  - Quadrant: —
  - GPU-amenable: maybe — pairwise alignment is SW-like
  - Upstream license: `BSD-2-Clause OR GPL-3`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-vsearch` with `--mode derep/chimera/cluster/pairalign` subcommands)
  - Consumes primitives: `noodles-fastq`, `noodles-fasta`, `rsomics-kmer`, `block-aligner`, `linfa-clustering`
  - Notes: Big rewrite — covers dereplication, chimera detection, pairwise alignment, clustering. Split into `--mode derep`, `--mode chimera`, `--mode pairalign` subcommands and tackle one at a time. SIMD-critical inner loops.

- [ ] **`swarm`** — single-linkage amplicon clustering without global cutoffs.
  - Reference impl: `C++` · [torognes/swarm](https://github.com/torognes/swarm) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — single-linkage neighbour expansion is SIMT-friendly
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-vsearch` (clustering mode)
  - Consumes primitives: `rsomics-kmer`, `petgraph`, `linfa-clustering`
  - Notes: Small focused codebase, clean Rust port target. Algorithm is elegant: iterative unweighted neighbor expansion. Good "starter" project for a new contributor.

- [ ] **`Deblur`** — sub-OTU denoising with positive read-error model.
  - Reference impl: `Python` · [biocore/deblur](https://github.com/biocore/deblur) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `DADA2`, `UNOISE3`
  - Parallelism: Python multiprocessing
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — error pattern subtraction is small
  - Upstream license: `BSD-3-Clause`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-dada2` (alternative denoiser mode)
  - Consumes primitives: `noodles-fastq`, `polars`
  - Notes: Less popular than DADA2 in 2026. Algorithm is straightforward (subtract error pattern from observed sequences). Port only if there's explicit user demand; DADA2 covers most needs.

- [ ] **`PICRUSt2`** — phylogeny-based functional prediction from 16S.
  - Reference impl: `Python` (wraps HMMER, EPA-NG, SEPP) · [picrust/picrust2](https://github.com/picrust/picrust2) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing + upstream tool pthreads
  - SIMD: inherits from HMMER / EPA-NG
  - Quadrant: —
  - GPU-amenable: maybe — HMM scoring SIMT-friendly
  - Upstream license: `GPL-3`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-picrust`)
  - Consumes primitives: future `rsomics-hmm`, future EPA-NG/`pplacer` substitute, `linfa` (HSP), `polars`
  - Notes: Phylogenetic placement (EPA-NG) + hidden-state prediction (HSP). HSP is small (~scikit-learn-style logistic regression); easy Rust port. Placement step inherits whatever solution we land on for GTDB-Tk/phylogeny.
