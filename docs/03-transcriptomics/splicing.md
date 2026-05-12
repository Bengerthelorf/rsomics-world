# Alternative splicing analysis

> Detection and quantification of differential splicing events from
> bulk RNA-seq.

## Scope

Two complementary approaches:

1. **Event-based** — enumerate canonical AS event types (SE, MXE, A3SS,
   A5SS, RI), compute PSI per event, test for ΔPSI between conditions.
   rMATS, SUPPA2, MAJIQ, Whippet, JUM.
2. **Junction-based** — model intron-junction usage as a multinomial,
   test for differential usage. LeafCutter.

Isoform-switch testing on the *transcript* level (Salmon abundances →
isoform fractions → tests) is also in scope (IsoformSwitchAnalyzeR).
Single-cell splicing is out of scope; see
[`../04-single-cell/`](../04-single-cell/).

## Design notes

- Splicing tools are mostly Python (rMATS, SUPPA2, MAJIQ) and R
  (LeafCutter, IsoformSwitchAnalyzeR). Whippet is Julia. There is no
  significant Rust presence yet.
- The performance bottleneck is junction-read counting from BAMs, which
  is a great match for `noodles-bam` + `petgraph` (splice graph) +
  `rayon`. A "fast junction matrix builder" Rust crate would be useful
  to *every* downstream tool here.
- Statistical kernels (rMATS hierarchical model, SUPPA2 empirical
  density, LeafCutter Dirichlet-multinomial) are small numerical
  routines that are reasonable Rust ports — but they need
  `ndarray-stats` + `statrs` to be in good shape first.
- The right abstraction for an rsomics splicing crate is a typed
  `SpliceGraph` over `petgraph`, an event enumerator, and a junction
  counter — feed any of the existing statistical methods from there.

## TODO

- [ ] **`rMATS`** (rMATS-turbo) — event-based ΔPSI from replicate
  RNA-seq.
  - Reference impl: `C++ / Python` · [Xinglab/rmats-turbo](https://github.com/Xinglab/rmats-turbo) · `BSD-3-Clause / mixed`
  - Existing Rust: none.
  - Existing non-C alternatives: SUPPA2 (Python, alignment-free).
  - Priority: `P0`
  - Notes: Most widely used event-based tool. A Rust rewrite would split
    naturally into (a) junction-counting (`noodles` + interval logic)
    and (b) the hierarchical-likelihood statistical kernel. The
    statistical model is in a published supplement and tractable.

- [ ] **`SUPPA2`** — alignment-free splicing analysis from transcript
  abundances.
  - Reference impl: `Python` · [comprna/SUPPA](https://github.com/comprna/SUPPA) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: none.
  - Priority: `P1`
  - Notes: Lightweight and pairs naturally with Salmon/kallisto/oarfish
    output. Small codebase, mostly NumPy; a clean Rust rewrite is
    feasible in a single crate (`rsomics-suppa`). Best splicing
    rewrite to attempt early.

- [ ] **`LeafCutter`** — intron-cluster differential splicing.
  - Reference impl: `C++ / R / Python` · [davidaknowles/leafcutter](https://github.com/davidaknowles/leafcutter) · `Apache-2.0`
  - Existing Rust: none.
  - Existing non-C alternatives: none.
  - Priority: `P1`
  - Notes: The intron-extraction step (`leafcutter_cluster_regtools.py`)
    is a clean Rust rewrite target — junction extraction from BAM,
    intron clustering. The Dirichlet-multinomial R kernel is small and
    can stay in R (extendr) or move to Rust via `statrs`.

- [ ] **`MAJIQ`** — local-splicing variations (LSV) detection and
  quantification.
  - Reference impl: `C++ / Python` · [biociphers/majiq](https://majiq.biociphers.org/) · `academic / non-commercial`
  - Existing Rust: none.
  - Existing non-C alternatives: none.
  - Priority: `P2`
  - Notes: License is restrictive (academic-only) so we **cannot** ship
    a derivative crate; clean-room reimplementation only. Even then,
    MAJIQ's algorithm is the most complex in this category and porting
    is a research-grade project.

- [~] **`Whippet`** — fast event-level splicing quantification using a
  splice graph.
  - Reference impl: `Julia` · [timbitz/Whippet.jl](https://github.com/timbitz/Whippet.jl) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: Whippet itself is Julia, the only
    non-C entry in this list.
  - Priority: `P2`
  - Notes: Marked `[~]` because the existing Julia implementation is
    already a clean modern rewrite of the same idea. A Rust port would
    primarily be motivated by deployability (no Julia runtime), not
    performance. Useful reference for splice-graph design in
    `rsomics-splice`.

- [ ] **`JUM`** — junction-usage modelling for splicing.
  - Reference impl: `Perl / R / Python` · [qqsong/JUM](https://github.com/qqsong/JUM) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: rMATS, LeafCutter.
  - Priority: `P2`
  - Notes: Niche tool. Listed for benchmark coverage; rewrite not
    planned.

- [ ] **`IsoformSwitchAnalyzeR`** — isoform-fraction-based switching
  analysis on top of Salmon / kallisto.
  - Reference impl: `R` · [Bioconductor IsoformSwitchAnalyzeR](https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html) · `GPL-2`
  - Existing Rust: none.
  - Existing non-C alternatives: none.
  - Priority: `P2`
  - Notes: Predominantly R wiring around DEXSeq / DRIMSeq, plus calls to
    domain-prediction tools. Right interop layer is `extendr`; no
    porting target on the rsomics side.
