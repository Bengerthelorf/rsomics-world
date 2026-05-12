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

- ML tree inference is dominated by **IQ-TREE2** and **RAxML-NG**. Both
  have rich model selection (ModelFinder / ModelTest-NG), bootstrap
  variants (UFBoot, SH-aLRT), partition models, and concordance-factor
  output. Reaching feature parity from Rust is multi-year. Phase 4+.
- Bayesian inference (MrBayes, BEAST2) is a research-software ecosystem
  in itself — partitioned models, prior families, MCMC samplers, calibrated
  dating. **BEAST2 is on the JVM (Java)** so we get nothing from rewriting
  in Rust except a faster MCMC; the model-specification DSL is the value.
- UShER is the right answer for ultra-large placements (SARS-CoV-2 scale).
  It's C++/SIMD; the algorithmic core (parsimony placement on a mutation-
  annotated tree) is published. Rust port is well-bounded but only useful
  if the upstream UShER stops being maintained.
- ASTRAL is Java. ASTER (the C++ rewrite, by Chao Zhang) is faster and
  more featureful; ASTRAL-Pro for paralogs lives in the same family.
  Pure-Rust ASTER would be a clean, small project.
- Tree IO (Newick, NEXUS, NeXML) is unsettled in Rust. `rust-newick`,
  `phylo-rs`, etc. exist but no clear winner. **Picking a canonical
  tree-data-structure crate is a Phase 1 task for this module.**
- License watch: RAxML-NG **AGPL-3** (note: this is more restrictive than
  GPL-3), IQ-TREE2 **GPL-2**, MrBayes **GPL-3**, BEAST2 **LGPL-2.1**,
  FastTree2 **GPL-2**, PhyML **CeCILL** (French GPL-compatible),
  ASTRAL **Apache-2.0**, ASTER **Apache-2.0** (check repo), UShER **MIT**
  (check), MEGA **proprietary-free**.

## TODO

- [ ] **`RAxML-NG`** — ML phylogeny inference (Next Generation rewrite of RAxML).
  - Reference impl: `C++` · [amkozlov/raxml-ng](https://github.com/amkozlov/raxml-ng) · `AGPL-3`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: AGPL-3 inheritance is the strictest in this module. Clean-room
    pure-Rust ML inference is the only path. RAxML-NG's strength is its
    pthreads-based parallel SPR search; Rust + `rayon` is well-suited.

- [ ] **`IQ-TREE2`** / `IQ-TREE3` — ML phylogenetics with ModelFinder and UFBoot.
  - Reference impl: `C++` · [iqtree/iqtree2](https://github.com/iqtree/iqtree2) · `GPL-2`
  - Existing Rust: none verified
  - Existing non-C alternatives: `RAxML-NG`
  - Priority: `P0`
  - Notes: The most-used ML inference tool in 2024-2026. Rust port is a
    huge undertaking but high-leverage. The model library (~50 substitution
    models + partition support + mixture models) is most of the work.
    Strategy: ship `rsomics-phylolikelihood` first (Felsenstein pruning +
    common DNA/AA models), then add tree-search heuristics on top.

- [ ] **`MrBayes`** — Bayesian phylogenetic inference (MCMC).
  - Reference impl: `C` · [NBISweden/MrBayes](https://github.com/NBISweden/MrBayes) · `GPL-3`
  - Existing Rust: none verified
  - Existing non-C alternatives: `BEAST2` (JVM, different ecosystem)
  - Priority: `P2`
  - Notes: Bayesian MCMC machinery is heavy; user demand has shifted
    toward BEAST2 and toward ML+UFBoot. Pure-Rust port is research-grade
    and probably not worth the effort vs. extending IQ-TREE-style approximations.

- [ ] **`BEAST` / `BEAST2`** — Bayesian Evolutionary Analysis by Sampling Trees.
  - Reference impl: `Java` · [CompEvol/beast2](https://github.com/CompEvol/beast2) · `LGPL-2.1`
  - Existing Rust: none verified
  - Existing non-C alternatives: `MrBayes`, `RevBayes`
  - Priority: `P2`
  - Notes: JVM-based, plugin-rich. Rust has no advantage over Java here
    on raw perf for MCMC. Skip.

- [ ] **`FastTree2`** — approximate ML for very large trees.
  - Reference impl: `C` · [morgannprice/fasttree](http://www.microbesonline.org/fasttree/) · `GPL-2`
  - Existing Rust: none verified
  - Existing non-C alternatives: `IQ-TREE2 --fast`
  - Priority: `P2`
  - Notes: Largely superseded by IQ-TREE2's fast mode and by UShER for
    ultra-large trees. Single C file (~16k LOC); a Rust port would be a
    weekend project but the use case is shrinking.

- [ ] **`MEGA`** — GUI phylogenetics for teaching/exploration.
  - Reference impl: `C++` (closed source for some versions) · [megasoftware.net](https://www.megasoftware.net/) · proprietary-free
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Closed source GUI. No port path; not a serious production
    target. Listed for completeness.

- [ ] **`PhyML`** — fast ML with NNI/SPR moves and aLRT support.
  - Reference impl: `C` · [stephaneguindon/phyml](https://github.com/stephaneguindon/phyml) · CeCILL (GPL-compatible)
  - Existing Rust: none verified
  - Existing non-C alternatives: `RAxML-NG`, `IQ-TREE2`
  - Priority: `P2`
  - Notes: Older alternative to RAxML-NG / IQ-TREE2. Skip; both successors
    cover this niche better.

- [x] **`UShER`** — ultra-fast sample placement onto existing trees.
  - Reference impl: `C++` · [yatisht/usher](https://github.com/yatisht/usher) · `MIT` (check repo)
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P0` (adopt; consider Rust port only if upstream becomes stale)
  - Notes: Built for SARS-CoV-2-scale (10M+-leaf) trees. Mutation-annotated
    tree (MAT) data structure is published; pure-Rust port is a clean,
    bounded project. Useful for both pathogen surveillance and general
    epi-phylogenetics. Current MIT licensing is friendly.

- [ ] **`TreeShrink`** — outlier-branch detection across tree sets.
  - Reference impl: `Python` (+ `R`) · [uym2/TreeShrink](https://github.com/uym2/TreeShrink) · check repo (GPL likely)
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Small Python codebase, well-defined algorithm. Small Rust port
    would slot neatly into `rsomics-phylo` as a tree-QC utility.

- [ ] **`ASTRAL`** / `ASTER` / `ASTRAL-Pro` — coalescent species-tree inference.
  - Reference impl: `Java` (ASTRAL) → `C++` (ASTER / ASTRAL-Pro) · [chaoszhang/ASTER](https://github.com/chaoszhang/ASTER) · check repo (Apache-2.0 reported)
  - Existing Rust: none verified
  - Existing non-C alternatives: `ASTER` is the C++ rewrite (preferred over Java ASTRAL)
  - Priority: `P1`
  - Notes: Target ASTER, not the Java ASTRAL. Algorithm (quartet score
    maximization on a constrained search space) is well-defined and a
    natural fit for `rayon`-parallel Rust. Medium-size project.
