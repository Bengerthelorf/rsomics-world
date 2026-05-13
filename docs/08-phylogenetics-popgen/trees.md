# Phylogenetic tree inference

> Maximum-likelihood, Bayesian, distance, coalescent species-tree, and
> ultra-large sample-placement tree inference and manipulation.

## Scope

Includes: ML tree inference (RAxML-NG, IQ-TREE 2/3, PhyML, FastTree2),
Bayesian inference (MrBayes, BEAST/BEAST2), heuristic/distance methods
(MEGA), coalescent species-tree (ASTRAL/ASTER family), tree placement at
scale (UShER), and tree-set quality control (TreeShrink). Excludes:
multiple-sequence alignment (see [alignment-msa](alignment-msa.md)) and
population-genetics analyses (see [population-genetics](population-genetics.md)).

## Design notes

- ML tree inference is dominated by **IQ-TREE2** and **RAxML-NG**.
- Bayesian inference (MrBayes, BEAST2) is a research-software ecosystem
  in itself. BEAST2 is JVM-based.
- UShER is the right answer for ultra-large placements (SARS-CoV-2 scale).
- ASTRAL is Java; ASTER (the C++ rewrite) is the modern target.
- Tree IO (Newick, NEXUS, NeXML) is unsettled in Rust — no clear canonical
  Rust crate yet. **Picking one is a Phase 1 task for this module.**
- License watch: RAxML-NG **AGPL-3** (strictest), IQ-TREE2 **GPL-2**,
  MrBayes **GPL-3**, BEAST2 **LGPL-2.1**, FastTree2 **GPL-2**, PhyML
  **CeCILL**, ASTRAL **Apache-2.0**, ASTER **Apache-2.0**, UShER **MIT**,
  MEGA **proprietary-free**.

## TODO

- [ ] **`RAxML-NG`** — ML phylogeny inference (Next Generation rewrite of RAxML).
  - Reference impl: `C++` · [amkozlov/raxml-ng](https://github.com/amkozlov/raxml-ng) · `AGPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads + MPI
  - SIMD: upstream SSE/AVX
  - Quadrant: —
  - GPU-amenable: maybe — phylogenetic likelihood is dense linear algebra
  - Upstream license: `AGPL-3`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-raxml`)
  - Consumes primitives: `noodles-fasta`, future `rsomics-phylo-tree` (tree data structure crate), future `rsomics-phylo-likelihood`, `ndarray-linalg`, `rayon`
  - Notes: AGPL-3 inheritance is the strictest in this module. Clean-room pure-Rust ML inference is the only path. Strength is pthreads-based parallel SPR search; Rust + `rayon` is well-suited.

- [ ] **`IQ-TREE2`** / `IQ-TREE3` — ML phylogenetics with ModelFinder and UFBoot.
  - Reference impl: `C++` · [iqtree/iqtree2](https://github.com/iqtree/iqtree2) · `GPL-2`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `RAxML-NG`
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE/AVX
  - Quadrant: —
  - GPU-amenable: maybe — Felsenstein pruning is dense linear algebra
  - Upstream license: `GPL-2`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-iqtree`)
  - Consumes primitives: future `rsomics-phylo-likelihood` (Felsenstein pruning), future `rsomics-phylo-tree`, `ndarray-linalg`, `rayon`
  - Notes: Most-used ML inference tool in 2024-2026. Rust port is a huge undertaking but high-leverage. Strategy: ship `rsomics-phylo-likelihood` first (Felsenstein pruning + common DNA/AA models), then add tree-search heuristics on top.

- [ ] **`MrBayes`** — Bayesian phylogenetic inference (MCMC).
  - Reference impl: `C` · [NBISweden/MrBayes](https://github.com/NBISweden/MrBayes) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `BEAST2` (JVM, different ecosystem)
  - Parallelism: upstream pthreads + MPI
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — MCMC parallelises across chains
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-iqtree` (Bayesian mode if added)
  - Consumes primitives: same as IQ-TREE2 plus MCMC machinery
  - Notes: Bayesian MCMC machinery is heavy; user demand has shifted toward BEAST2 and ML+UFBoot. Pure-Rust port is research-grade.

- [ ] **`BEAST` / `BEAST2`** — Bayesian Evolutionary Analysis by Sampling Trees.
  - Reference impl: `Java` · [CompEvol/beast2](https://github.com/CompEvol/beast2) · `LGPL-2.1`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `MrBayes`, `RevBayes`
  - Parallelism: JVM threading
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — MCMC parallelisation
  - Upstream license: `LGPL-2.1`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: JVM-based, plugin-rich. Rust has no advantage over Java here on raw perf for MCMC. Skip.

- [ ] **`FastTree2`** — approximate ML for very large trees.
  - Reference impl: `C` · [microbesonline.org/fasttree](http://www.microbesonline.org/fasttree/) · `GPL-2`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `IQ-TREE2 --fast`
  - Parallelism: upstream OpenMP
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — approximate ML scoring
  - Upstream license: `GPL-2`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-iqtree` (fast mode)
  - Consumes primitives: same as IQ-TREE2
  - Notes: Largely superseded by IQ-TREE2's fast mode and by UShER for ultra-large trees. Single C file (~16k LOC); a Rust port would be a weekend project but the use case is shrinking.

- [ ] **`MEGA`** — GUI phylogenetics for teaching/exploration.
  - Reference impl: `C++` · [megasoftware.net](https://www.megasoftware.net/) · proprietary-free
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: closed
  - SIMD: closed
  - Quadrant: —
  - GPU-amenable: no — GUI
  - Upstream license: proprietary-free
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Closed source GUI. No port path; not a serious production target.

- [ ] **`PhyML`** — fast ML with NNI/SPR moves and aLRT support.
  - Reference impl: `C` · [stephaneguindon/phyml](https://github.com/stephaneguindon/phyml) · `CeCILL` (GPL-compatible)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `RAxML-NG`, `IQ-TREE2`
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — same as RAxML
  - Upstream license: `CeCILL`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Older alternative to RAxML-NG / IQ-TREE2. Skip; both successors cover this niche better.

- [A] **`UShER`** — ultra-fast sample placement onto existing trees.
  - Reference impl: `C++` · [yatisht/usher](https://github.com/yatisht/usher) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — parsimony placement parallelises trivially
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `subcommand-of-rsomics-iqtree` (wraps the upstream UShER binary via subprocess); a Rust port to `rsomics-usher` is a fallback if upstream becomes stale
  - Consumes primitives: subprocess call to UShER binary; `noodles-vcf` for input formatting
  - Notes: Built for SARS-CoV-2-scale (10M+-leaf) trees. MIT, actively maintained. The `[A]` mark per the post-Phase-1 schema revision: adopt the upstream binary via subprocess, no rewrite planned. The MAT data structure is published; pure-Rust port remains a clean fallback option if upstream stalls.

- [ ] **`TreeShrink`** — outlier-branch detection across tree sets.
  - Reference impl: `Python` (+ `R`) · [uym2/TreeShrink](https://github.com/uym2/TreeShrink) · `GPL`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: no — small algorithm
  - Upstream license: `GPL`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-iqtree` (tree-QC mode)
  - Consumes primitives: future `rsomics-phylo-tree`
  - Notes: Small Python codebase, well-defined algorithm. Small Rust port would slot neatly into `rsomics-phylo` as a tree-QC utility.

- [ ] **`ASTRAL`** / `ASTER` / `ASTRAL-Pro` — coalescent species-tree inference.
  - Reference impl: `Java` (ASTRAL) → `C++` (ASTER / ASTRAL-Pro) · [chaoszhang/ASTER](https://github.com/chaoszhang/ASTER) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `ASTER` is the C++ rewrite (preferred over Java ASTRAL)
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — quartet score evaluation parallelises
  - Upstream license: `Apache-2.0`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-aster`)
  - Consumes primitives: future `rsomics-phylo-tree`, `rayon`
  - Notes: Target ASTER, not the Java ASTRAL. Algorithm (quartet score maximization on a constrained search space) is well-defined and a natural fit for `rayon`-parallel Rust.
