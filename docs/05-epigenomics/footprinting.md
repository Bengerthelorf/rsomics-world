# Transcription-factor footprinting

> Inferring TF binding from ATAC-seq / DNase-seq signal at base-pair
> resolution.

## Scope

De novo and motif-guided footprinting tools: TOBIAS, HINT-ATAC, BinDNase,
PIQ, CENTIPEDE, plus the V-plot family of visualisation / characterisation
tools.

Upstream peak calling and signal-track generation are in
[`peak-calling.md`](peak-calling.md). Motif-database / scanning crates
(`MEME`, `FIMO`, etc.) are not in scope here — they belong with motif
analysis under a future regulatory-genomics topic.

## Design notes

- Footprinting is a small, mathematically tight sub-area: per-position
  Tn5 / DNase cut counts (a coverage track) + a bias model + a per-
  motif scoring. The data pipeline (BAM → bias-corrected per-base
  coverage) is exactly the kind of IO-heavy work where Rust + `noodles`
  is a natural win.
- TOBIAS is the current state-of-the-art and is pure Python + Cython;
  HINT-ATAC is older Python; PIQ and CENTIPEDE are R.
- An `rsomics-footprint` crate that provides (a) Tn5 bias correction,
  (b) per-motif footprint scoring, (c) differential binding statistics
  would consolidate this area. The bias-correction model is the
  largest piece; the rest is small.
- This is a Phase-3 niche — high value to chromatin researchers, lower
  total user count than peak calling.

## TODO

- [ ] **`TOBIAS`** — Tn5-bias-corrected ATAC-seq footprinting framework.
  - Reference impl: `Python / Cython` · [loosolab/TOBIAS](https://github.com/loosolab/TOBIAS) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — per-position bias correction is SIMT-friendly
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-footprint`)
  - Consumes primitives: `noodles-bam`, `ndarray`, `rsomics-coverage`, future `rsomics-motif`, `rsomics-intervals`
  - Notes: Algorithm is well-documented (ATACorrect bias model + BINDetect differential binding). Compact Rust rewrite using `noodles-bam`, `ndarray`, and a motif-scanning primitive. Best Rust target in this sub-area.

- [ ] **`HINT-ATAC`** — HMM-based footprinting for ATAC-seq.
  - Reference impl: `Python` · [CostaLab/reg-gen](https://github.com/CostaLab/reg-gen) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: maybe — HMM Baum-Welch is GPU-friendly at scale
  - Upstream license: `GPL-3.0`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-footprint`
  - Consumes primitives: `noodles-bam`, `rsomics-coverage`, future `rsomics-stats` (HMM machinery)
  - Notes: Older HMM footprinter; TOBIAS outperforms it in benchmarks. Listed for completeness.

- [ ] **`BinDNase`** — DNase-seq footprinter using bin-level statistics.
  - Reference impl: `Python / R` · academic publication · unspecified
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: limited
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — niche, no upside
  - Upstream license: unspecified
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Niche; rarely used today.

- [ ] **`PIQ`** — Bayesian footprinter (Pique).
  - Reference impl: `R` · [thashim/piq-single](https://bitbucket.org/thashim/piq-single/) · unspecified (Bitbucket-hosted, GitHub aliveness check N/A)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: R BiocParallel
  - SIMD: BLAS
  - Quadrant: —
  - GPU-amenable: no — niche, older
  - Upstream license: unspecified
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Niche; older method. Not a porting target. Hosted on Bitbucket; aliveness check via gh not possible.

- [ ] **`CENTIPEDE`** — original Bayesian footprinter.
  - Reference impl: `R` · CRAN · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: R BiocParallel
  - SIMD: BLAS
  - Quadrant: —
  - GPU-amenable: no — classical
  - Upstream license: `GPL-3.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Classical method; superseded by TOBIAS in practice.

- [ ] **ATAC V-plot tooling** — fragment-size × position visualisation tools (Greenleaf-lab style `ATAC-Vsignal`, `pyatac`).
  - Reference impl: `Python` · canonical Greenleaf-lab `pyatac` URL is dead (no live GitHub repo found via search). V-plot concept is widely used; implementations scattered across lab repos · MIT-style typical
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: trivial parallel across features
  - SIMD: auto-vectorize
  - Quadrant: —
  - GPU-amenable: maybe — 2D histogram of fragment midpoint × length is GPU-trivial
  - Upstream license: MIT (typical for ATAC analysis tooling)
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-footprint` (V-plot mode)
  - Consumes primitives: `noodles-bam`, `ndarray`, `rsomics-intervals`
  - Notes: Small pure-Rust win — the V-plot is just a 2D histogram of fragment midpoint × fragment length around a feature, parallelised across features. Bundle into `rsomics-footprint` as a side feature. **Original entry's `GreenleafLab/pyatac` URL is dead** — logged to `.autopilot/needs-review/external-2026-05-14.md`; will need user adjudication on whether to find a successor V-plot tool or treat the concept as standalone.
