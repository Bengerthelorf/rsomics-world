# ChIP-seq / ATAC-seq pipelines and chromatin state

> Full ChIP-seq / ATAC-seq pipelines (QC, trimming, alignment, peak
> calling, signal tracks) and downstream chromatin-state segmentation.

## Scope

End-to-end pipelines (nf-core/atacseq, nf-core/chipseq, ENCODE ATAC),
specialised peak callers shipped only within those pipelines (Genrich),
and chromatin-state segmentation (ChromHMM, Segway).

Standalone peak callers live in [`peak-calling.md`](peak-calling.md).
TF footprinting is in [`footprinting.md`](footprinting.md).

## Design notes

- Real-world ChIP / ATAC analysis is orchestrated — every step is a
  separate binary glued together by Snakemake (encode-dcc, AIAP) or
  Nextflow (nf-core). The rsomics deliverable here is *less about
  rewriting algorithms* and more about producing high-quality Rust
  binaries (trimmer, aligner, dedup, peak caller, signal tracks) that
  these pipelines can adopt.
- ChromHMM and Segway are unusual entries in this module: they are
  *segmentation* tools that consume per-bin ChIP signal across multiple
  marks and emit a discrete annotation. Both have small, mathematically
  clean cores (multivariate HMM / DBN) and are good Rust targets.
- Genrich is a pure-C peak caller bundled inside ATAC pipelines (its
  ATAC-mode is the default in nf-core/atacseq replicates). Rust-port
  priority is low only because MACS3 with ATAC-mode covers most use
  cases.

## TODO

- [ ] **`ENCODE ATAC-seq pipeline`** — reference ATAC-seq pipeline.
  - Reference impl: `WDL / Python / Shell` · [ENCODE-DCC/atac-seq-pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: nf-core/atacseq (Nextflow)
  - Parallelism: WDL workflow engine (Cromwell)
  - SIMD: inherited from invoked binaries
  - Quadrant: —
  - GPU-amenable: no — orchestration layer
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: —
  - Consumes primitives: future rsomics-* binaries
  - Notes: Orchestration only — rewriting the pipeline is out of scope for rsomics. The contribution we make is providing Rust components (trimmer, aligner, MACS-rs) that ENCODE can adopt.

- [ ] **`nf-core/atacseq`** — community Nextflow ATAC-seq pipeline.
  - Reference impl: `Nextflow / Python` · [nf-core/atacseq](https://github.com/nf-core/atacseq) · `MIT`
  - Existing Rust: none of the pipeline itself; many tools it invokes have potential Rust replacements
  - Existing Rust kind: `none`
  - Existing non-C alternatives: ENCODE ATAC pipeline
  - Parallelism: Nextflow workflow engine
  - SIMD: inherited from invoked binaries
  - Quadrant: —
  - GPU-amenable: no — orchestration layer
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: —
  - Consumes primitives: future rsomics-* binaries (rsomics-macs, rsomics-bwa, rsomics-fastp etc.)
  - Notes: Track which Rust binaries can be swapped in. Pipeline orchestration is out of scope here; see module 09.

- [ ] **`nf-core/chipseq`** — community Nextflow ChIP-seq pipeline.
  - Reference impl: `Nextflow / Python` · [nf-core/chipseq](https://github.com/nf-core/chipseq) · `MIT`
  - Existing Rust: as above
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Nextflow workflow engine
  - SIMD: inherited
  - Quadrant: —
  - GPU-amenable: no — orchestration layer
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: —
  - Consumes primitives: future rsomics-* binaries
  - Notes: Same notes as `nf-core/atacseq`.

- [ ] **`AIAP`** — ATAC-seq Integrative Analysis Pipeline.
  - Reference impl: `Python / Shell` · [Zhang-lab/ATAC-seq_QC_analysis](https://github.com/Zhang-lab/ATAC-seq_QC_analysis) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: shell orchestration
  - SIMD: inherited
  - Quadrant: —
  - GPU-amenable: no — orchestration layer
  - Upstream license: `MIT`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Smaller user base than nf-core/atacseq. Listed for completeness.

- [ ] **`Genrich`** — pure-C peak caller bundled in ATAC pipelines.
  - Reference impl: `C` · [jsh58/Genrich](https://github.com/jsh58/Genrich) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: MACS3 ATAC mode
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — pileup + replicate p-value combination
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-macs` (ATAC replicate-aware mode)
  - Consumes primitives: `noodles-bam`, `rsomics-coverage`, `statrs`, `rsomics-intervals`
  - Notes: Cleanly written C, ~5 kLoC, replicate-aware p-value combination is its distinguishing feature. Small, clean Rust port target; reuse the `rsomics-coverage` primitive.

- [ ] **`ChromHMM`** — multivariate HMM chromatin-state segmentation.
  - Reference impl: `Java` · [jernst98/ChromHMM](https://github.com/jernst98/ChromHMM) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: ChromHMM itself (Java)
  - Parallelism: JVM threading
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: maybe — Baum-Welch can run on GPU at scale
  - Upstream license: `GPL-3.0`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-chromhmm`)
  - Consumes primitives: `ndarray-stats`, `statrs`, future `rsomics-stats` (HMM machinery)
  - Notes: Pure HMM with Bernoulli emissions on binarised marks. Small algorithmic core (Baum-Welch). `ndarray-stats` covers the math. A Rust ChromHMM would be a clean focused crate and could share the HMM machinery with future short-read aligners (HISAT-style splice models). Reference URL in the original entry was wrong (`ernstlab/ChromHMM` → 404); corrected to `jernst98/ChromHMM` (Jason Ernst's personal repo) after a GitHub search.

- [ ] **`Segway`** — DBN-based chromatin-state segmentation.
  - Reference impl: `Python / C++ (GMTK)` · [hoffmangroup/segway](https://github.com/hoffmangroup/segway) · `GPL-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: ChromHMM (less expressive but simpler)
  - Parallelism: upstream Python + GMTK
  - SIMD: GMTK's hand SIMD
  - Quadrant: —
  - GPU-amenable: maybe — DBN inference has GPU variants
  - Upstream license: `GPL-2.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Dependency on the GMTK dynamic Bayesian network library makes a Rust port large. Lower priority than ChromHMM.
