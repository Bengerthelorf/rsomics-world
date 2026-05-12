# 05 — Epigenomics

> Peak calling, chromatin pipelines, DNA methylation, 3D chromatin
> conformation, and transcription-factor footprinting.

## Sub-topics

- [`peak-calling.md`](peak-calling.md) — MACS2/3, HOMER findPeaks, SEACR,
  SICER, PeakSeq, GoPeaks, GEM.
- [`chip-atac-pipelines.md`](chip-atac-pipelines.md) — ENCODE ATAC-seq
  pipeline, nf-core/atacseq, nf-core/chipseq, AIAP, Genrich, ChromHMM,
  Segway.
- [`methylation.md`](methylation.md) — Bismark, BWA-meth, methylKit,
  methylpy, MethylDackel, modkit, nanopolish.
- [`chromatin-3d.md`](chromatin-3d.md) — cooler, cooltools, HiC-Pro,
  HiCExplorer, Juicer, pairtools, distiller, FAN-C, chromosight,
  mustache.
- [`footprinting.md`](footprinting.md) — TOBIAS, HINT-ATAC, BinDNase,
  PIQ, CENTIPEDE, ATAC-V-plot tools.

## Cross-cutting design notes

- Bulk epigenomics is one of the most "C-and-Python" sub-fields of
  bioinformatics. Almost every tool here is a C / C++ binary with a
  Python or R wrapper — fertile ground for Rust ports.
- The single Rust success story is **`modkit`** (Oxford Nanopore), the
  reference tool for modified-base summarisation on long reads. Adopt
  as-is.
- Single-cell scATAC analysis is covered in
  [`../04-single-cell/`](../04-single-cell/), particularly via
  `SnapATAC2` (Rust + Python). The peak-calling / fragment-handling
  primitives here should be shared.
- 3D chromatin (Hi-C, Micro-C, HiChIP) is dominated by Python + HDF5
  (cooler) and Java (Juicer). The pairtools + cooler stack is a clean
  Rust port target because the data model (`pairs` text format, `cool`
  HDF5) is small and well-specified.
- DNA methylation has two regimes:
  - **Short-read bisulfite** (Bismark, BWA-meth) — algorithmically
    similar to short-read alignment, so a Rust solution rides on the
    short-read aligner work from module 02.
  - **Long-read modified-base calling** (modkit, nanopolish) — already
    mostly Rust via modkit; only the legacy nanopolish remains C++.
