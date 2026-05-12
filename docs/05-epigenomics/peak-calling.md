# Peak calling

> Identification of enriched genomic regions from ChIP-seq, ATAC-seq,
> CUT&RUN, and CUT&Tag data.

## Scope

Algorithms that take aligned reads (BAM / BED) and emit a list of
enriched intervals plus statistics. Downstream integration with
chromatin-state segmentation (ChromHMM, Segway) and full pipelines
(nf-core/chipseq) is in [`chip-atac-pipelines.md`](chip-atac-pipelines.md).
Hi-C "peak"/loop callers (chromosight, mustache, HiCCUPS) live in
[`chromatin-3d.md`](chromatin-3d.md).

## Design notes

- MACS2/3 is essentially universal — the rsomics goal is a Rust binary
  whose output BED / narrowPeak / broadPeak is bit-identical to
  `macs3 callpeak` for the typical defaults so downstream pipelines
  (nf-core/chipseq, ENCODE) work unchanged.
- The bottleneck in MACS is the genome-wide Poisson background
  estimation and the sliding-window peak shape model — both naturally
  parallelise across chromosomes with `rayon`, and `noodles-bam` reads
  the input cheaply.
- Several peak callers are already in modern languages: `GoPeaks` is
  Go, `GEM` is Java, `epic2` is Cython. Rust-rewrite priority should
  favour tools that are pure-C and pure-Python (MACS, HOMER, SEACR,
  SICER) where Rust has the largest performance / safety win.
- All peak callers share the same upstream primitive — a coverage
  pileup from BAM with optional fragment-length extension or paired-end
  fragment reconstruction. `rsomics-coverage` is the right shared crate.

## TODO

- [ ] **`MACS2`** — model-based ChIP-seq peak caller (legacy).
  - Reference impl: `Python / C` · [macs3-project/MACS](https://github.com/macs3-project/MACS) · `BSD-3-Clause`
  - Existing Rust: none.
  - Existing non-C alternatives: GoPeaks (Go), epic2 (Cython).
  - Priority: `P1`
  - Notes: Most pipelines now pin MACS3. Keep MACS2 listed for
    reproducibility of older studies; the rsomics target is MACS3 with
    a compatibility flag.

- [ ] **`MACS3`** — current MACS with ATAC-seq, scATAC, and broad-peak
  improvements.
  - Reference impl: `Python / Cython` · [macs3-project/MACS](https://github.com/macs3-project/MACS) · `BSD-3-Clause`
  - Existing Rust: none.
  - Existing non-C alternatives: GoPeaks (Go), epic2 (Cython).
  - Priority: `P0`
  - Notes: **Highest-impact Rust target in this module.** Pure-Rust
    `rsomics-macs` should match `narrowPeak` / `broadPeak` /
    `xls` outputs byte-for-byte at default settings. `noodles-bam` +
    `rsomics-coverage` + `statrs` (Poisson, log-gamma) covers it.

- [ ] **`HOMER findPeaks`** — tag-directory-based peak caller with
  GRO-seq / 4C / etc. modes.
  - Reference impl: `Perl / C++` · [HOMER](http://homer.ucsd.edu/homer/) · `unspecified open source`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: HOMER's tag directory is a non-standard format that downstream
    HOMER tools assume — a Rust port should preserve the format, not
    just the peak BEDs. Lower priority than MACS3 but still widely
    used.

- [ ] **`SEACR`** — Sparse Enrichment Analysis for CUT&RUN.
  - Reference impl: `R / Bash` · [FredHutch/SEACR](https://github.com/FredHutch/SEACR) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Small algorithm — global-distribution-based threshold over
    contiguous blocks. Compact Rust crate, < 1k LoC. Important for
    CUT&RUN pipelines.

- [ ] **`SICER` / `SICER2`** — spatial-clustering peak caller for
  broad histone marks.
  - Reference impl: `Python / C` · [zanglab/SICER2](https://github.com/zanglab/SICER2) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: `epic2` (Cython, much faster).
  - Priority: `P1`
  - Notes: epic2 is already a 30× speedup reimplementation. A Rust port
    of epic2's algorithm would be a small additional win and slot into
    `rsomics-peaks` cleanly.

- [ ] **`PeakSeq`** — original ENCODE-era ChIP-seq peak caller.
  - Reference impl: `C / Perl` · [gersteinlab/PeakSeq](https://github.com/gersteinlab/PeakSeq) · `GPL-2`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Largely historical. Listed for ENCODE legacy reproducibility.
    No porting target.

- [~] **`GoPeaks`** — modern Go peak caller for CUT&Tag.
  - Reference impl: `Go` · [maxsonBraunLab/gopeaks](https://github.com/maxsonBraunLab/gopeaks) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: GoPeaks itself (Go).
  - Priority: `P2`
  - Notes: Marked `[~]` because the existing Go implementation is
    already a clean modern rewrite. Listed for completeness; a Rust
    port is a low-priority duplication.

- [ ] **`GEM`** — Java peak caller with motif-aware refinement.
  - Reference impl: `Java` · [gem](http://groups.csail.mit.edu/cgs/gem/) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: GEM itself (Java).
  - Priority: `P2`
  - Notes: Niche; motif-coupled peak refinement is its differentiator.
    Not a porting target.
