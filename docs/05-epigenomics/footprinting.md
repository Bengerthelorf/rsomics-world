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
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Algorithm is well-documented (ATACorrect bias model + BINDetect
    differential binding). Compact Rust rewrite using `noodles-bam`,
    `ndarray`, and a motif-scanning primitive. Best Rust target in
    this sub-area.

- [ ] **`HINT-ATAC`** — HMM-based footprinting for ATAC-seq.
  - Reference impl: `Python` · [CostaLab/reg-gen](https://github.com/CostaLab/reg-gen) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Older HMM footprinter; TOBIAS outperforms it in benchmarks.
    Listed for completeness.

- [ ] **`BinDNase`** — DNase-seq footprinter using bin-level statistics.
  - Reference impl: `Python / R` · academic publication · `unspecified`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Niche; rarely used today.

- [ ] **`PIQ`** — Bayesian footprinter (Pique).
  - Reference impl: `R` · [vplaboratory/piq](https://bitbucket.org/thashim/piq-single/) · `unspecified`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Niche; older method. Not a porting target.

- [ ] **`CENTIPEDE`** — original Bayesian footprinter.
  - Reference impl: `R` · [Pique CRAN](https://cran.r-project.org/) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Classical method; superseded by TOBIAS in practice.

- [ ] **ATAC V-plot tooling** — fragment-size × position visualisation
  tools (Greenleaf-lab style `ATAC-Vsignal`, `pyatac`).
  - Reference impl: `Python` · [GreenleafLab/pyatac](https://github.com/GreenleafLab/pyatac) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Small pure-Rust win — the V-plot is just a 2D histogram of
    fragment midpoint × fragment length around a feature, parallelised
    across features. Bundle into `rsomics-footprint` as a side feature.
