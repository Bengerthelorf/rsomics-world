# Differential expression analysis

> Statistical testing of count / abundance changes between conditions for
> bulk RNA-seq.

## Scope

Methods that take a count matrix (or transcript abundance + bootstraps,
in the case of sleuth) and a sample-condition design matrix and emit a
table of per-gene / per-transcript log fold changes, statistics, and
adjusted p-values.

Single-cell DE belongs in [`../04-single-cell/analysis-core.md`](../04-single-cell/analysis-core.md).
Splicing-level DE is in [`splicing.md`](splicing.md). The actual matrix
construction is in [`quantification.md`](quantification.md).

## Design notes

- This area is firmly R / Bioconductor (DESeq2, edgeR, limma) plus a
  handful of Python ports (pyDESeq2). The methods are decades of careful
  empirical-Bayes / GLM work; reimplementing them is a research project,
  not a port.
- **Interop first.** Any rsomics-DE crate should plumb a Rust count
  matrix → `extendr` / `rextendr` bridge → DESeq2 / edgeR in R, with
  optional PyO3 path to pyDESeq2. Polars or Arrow IPC for the matrix
  layer keeps the boundary cheap.
- The kernels that are bottlenecks in pure-R / pure-Python — IRLS for the
  GLM, dispersion estimation, voom precision-weight calculation — are
  amenable to rewriting in Rust with `ndarray-linalg` + `rayon`. Wrap
  them behind the R / Python facade so users don't switch packages.
- `polars` is the right home for the result tables (gene × stats) and
  for the design-matrix construction utilities.
- Tools that ship their own quantifier upstream (sleuth + kallisto) need
  the bootstrap-aware variance model — that is genuinely novel work, not
  a port of an existing Rust kernel.

## TODO

- [ ] **`DESeq2`** — negative binomial GLM with shrinkage, the most-cited bulk RNA-seq DE method.
  - Reference impl: `R / C++` · [Bioconductor DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) · `LGPL-3+`
  - Existing Rust: none verified. Python port `pydeseq2` exists from owkin (not Rust)
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `pydeseq2` (Python)
  - Parallelism: R BiocParallel + C++ inner loops
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — IRLS for the GLM is dense linear algebra and SIMT-friendly
  - Upstream license: `LGPL-3+`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-de` as the umbrella; DESeq2 path is an interop bridge first, kernel-level rewrite later)
  - Consumes primitives: `polars`, `extendr`-bridge, future `rsomics-stats` (IRLS + dispersion), `ndarray-linalg`
  - Notes: A from-scratch Rust port is a research project. Realistic rsomics deliverable: a `rsomics-de` crate that produces DESeq2's input format (counts + colData), invokes DESeq2 via `extendr`, returns results as `polars` frames. Kernel-level rewrite of the `nbinomWaldTest` IRLS is a Phase-3 stretch goal.

- [ ] **`edgeR`** — empirical Bayes negative binomial GLM, original bulk RNA-seq DE method.
  - Reference impl: `R / C++` · [Bioconductor edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) · `GPL-2+`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: none widely-used
  - Parallelism: R BiocParallel
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — same dense-linear-algebra rationale as DESeq2
  - Upstream license: `GPL-2+`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-de` (one umbrella binary with `--method deseq2` / `--method edger` / etc.)
  - Consumes primitives: `polars`, `extendr`-bridge, future `rsomics-stats`, `ndarray-linalg`
  - Notes: Same interop strategy as DESeq2. `edgeR-v4`'s quasi-likelihood F-test is the current default; ensure any wrapper exposes it.

- [ ] **`limma-voom`** — precision-weighted linear models for RNA-seq.
  - Reference impl: `R` · [Bioconductor limma](https://bioconductor.org/packages/release/bioc/html/limma.html) · `GPL-2+`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: R BiocParallel
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: maybe — voom weights are per-gene independent
  - Upstream license: `GPL-2+`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-de`
  - Consumes primitives: `polars`, `extendr`-bridge, future `rsomics-stats`, `ndarray-linalg`
  - Notes: limma's empirical Bayes moderation is decades of careful statistics. The `voom` weight calculation itself is small and a plausible Rust kernel target; the moderated t-test machinery is not.

- [ ] **`sleuth`** — bootstrap-aware DE for kallisto / Salmon abundances.
  - Reference impl: `R` · [pachterlab/sleuth](https://github.com/pachterlab/sleuth) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `swish` (in `fishpond`) is the spiritual successor for Salmon
  - Parallelism: R BiocParallel
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: maybe — bootstrap variance decomposition parallelises trivially
  - Upstream license: `GPL-3.0`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-de` (a `--bootstrap-aware` mode)
  - Consumes primitives: `polars`, `extendr`-bridge, future `rsomics-stats`, `rsomics-kallisto` (for the bootstrap-emitting quantifier upstream)
  - Notes: Pairs naturally with a Rust kallisto port (see `quantification.md`). The technical-variance / biological-variance decomposition is the novel piece; would be a clean Rust crate if we have a bootstrap-emitting quantifier upstream.

- [ ] **`NOISeq`** — non-parametric DE for low-replicate experiments.
  - Reference impl: `R` · [Bioconductor NOISeq](https://www.bioconductor.org/packages/release/bioc/html/NOISeq.html) · `Artistic-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: R BiocParallel
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — ranking + percentile work
  - Upstream license: `Artistic-2.0`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-de`
  - Consumes primitives: `polars`, `extendr`-bridge
  - Notes: Niche. Pure-R, small codebase. Could be rewritten in Rust via `polars` + ranking utilities but low pipeline impact.

- [ ] **`EBSeq`** — empirical-Bayes hierarchical model for isoform DE.
  - Reference impl: `R` · [Bioconductor EBSeq](https://bioconductor.org/packages/release/bioc/html/EBSeq.html) · `Artistic-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: R BiocParallel
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: maybe — hierarchical EM is dense
  - Upstream license: `Artistic-2.0`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-de`
  - Consumes primitives: `polars`, `extendr`-bridge, RSEM-style input
  - Notes: Often paired with RSEM output. Niche. Wrapping via `extendr` is sufficient.
