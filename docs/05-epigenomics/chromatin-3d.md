# 3D chromatin (Hi-C, Micro-C, HiChIP)

> From read pairs to contact matrices, normalised maps, and topological
> feature calls (TADs, loops, compartments).

## Scope

The full Hi-C analysis stack:

1. **Read pair extraction** — pairtools / distiller.
2. **Matrix storage and manipulation** — cooler (.cool / .mcool),
   Juicer (.hic).
3. **Pipelines** — HiC-Pro, distiller, HiCExplorer, Juicer.
4. **Feature calling** — loops (chromosight, mustache, hicFindPeaks),
   TADs / compartments (cooltools, HiCExplorer), differential analysis
   (FAN-C / fanc-tools).

scHi-C lives in module 04 (spatial / single-cell).

## Design notes

- The Hi-C stack is unusually well-specified: the `pairs` text format
  is a 4DN standard; the `.cool` format is HDF5 with a tiny schema;
  the `.hic` format is open. Rust IO crates for these are missing but
  small — a focused `cooler-rs` would unlock the rest.
- The bottleneck in real-world Hi-C is the **read-pair extraction**
  step (pairtools parse + dedup) on hundreds of millions of read
  pairs. Pure-Rust pairtools — `noodles-bam` reader + dedup hash +
  parallel sort — could be 5–10× faster than the Python original.
- Loop / TAD callers are mostly small image-processing-on-Hi-C-matrix
  algorithms (chromosight uses pattern matching; mustache uses
  scale-space blob detection). `ndarray` + `imageproc` style Rust
  crates cover the math.
- Juicer (Java) is fast but closed-feeling; the open ecosystem (cooler
  + cooltools) is what the field is converging on. Rust effort should
  follow that trajectory.

## TODO

- [ ] **`cooler`** — HDF5-backed Hi-C contact matrix format and library.
  - Reference impl: `Python` · [open2c/cooler](https://github.com/open2c/cooler) · `BSD-3-Clause`
  - Existing Rust: none. `hdf5-rust` provides the IO primitives.
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Build `cooler-rs` first — every other 3D chromatin entry
    depends on reading/writing `.cool` / `.mcool`. Schema is tiny
    (a few HDF5 groups). One-week project for a focused Rust dev.

- [ ] **`cooltools`** — analysis on top of cooler (eigenvector
  compartments, insulation, dot calling).
  - Reference impl: `Python` · [open2c/cooltools](https://github.com/open2c/cooltools) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Once `cooler-rs` exists, each analytic in cooltools
    (insulation, eigendecomposition for A/B compartments, expected
    decay) is a small `ndarray` routine. Build piecewise.

- [ ] **`HiC-Pro`** — established Hi-C processing pipeline (mapping +
  filtering + matrix construction).
  - Reference impl: `Python / Bash` · [nservant/HiC-Pro](https://github.com/nservant/HiC-Pro) · `BSD-3-Clause`
  - Existing Rust: none.
  - Existing non-C alternatives: distiller (Nextflow).
  - Priority: `P1`
  - Notes: The "official" pipeline orchestration here is out of scope;
    the rsomics contribution is Rust components (rsomics-pairtools,
    rsomics-cooler) it can call.

- [ ] **`HiCExplorer`** — Hi-C visualisation and analysis suite.
  - Reference impl: `Python` · [deeptools/HiCExplorer](https://github.com/deeptools/HiCExplorer) · `GPL-3.0`
  - Existing Rust: none.
  - Priority: `P2`
  - Notes: User-facing tool; not a priority for Rust rewrite.

- [ ] **`Juicer` / `JuicerTools`** — Aiden lab pipeline + `.hic` toolkit.
  - Reference impl: `Java + C++` · [aidenlab/juicer](https://github.com/aidenlab/juicer) · `MIT`
  - Existing Rust: none. `hicstraw`-style readers exist but only in
    Python / C++.
  - Existing non-C alternatives: cooler / cooltools is the open
    counterpart.
  - Priority: `P2`
  - Notes: A Rust `.hic` reader/writer would be useful for legacy
    compatibility, but the field is converging on `.cool`.

- [ ] **`pairtools`** — extract `.pairs` from SAM/BAM, dedup, filter.
  - Reference impl: `Python / Cython` · [open2c/pairtools](https://github.com/open2c/pairtools) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: **High-impact Rust target.** Read-pair extraction is the
    pipeline bottleneck. `noodles-bam` reader + parallel dedup
    (xxhash3 → sharded hashmap) + parallel sort → `.pairs.gz` output.
    Match pairtools' format exactly so downstream cooler can read it.

- [ ] **`distiller`** — Nextflow Hi-C pipeline maintained by Open2C.
  - Reference impl: `Nextflow / Python` · [open2c/distiller-nf](https://github.com/open2c/distiller-nf) · `MIT`
  - Existing Rust: none.
  - Priority: `P2`
  - Notes: Pipeline orchestration; out of scope. Components above
    serve as drop-in replacements.

- [ ] **`hicFindPeaks`** (HOMER) — Hi-C loop caller within HOMER.
  - Reference impl: `Perl / C++` · part of HOMER · `unspecified`
  - Existing Rust: none.
  - Priority: `P2`
  - Notes: Niche; users now prefer chromosight or mustache. Listed for
    completeness.

- [ ] **`FAN-C` / `fanc-tools`** — multi-resolution Hi-C analysis suite.
  - Reference impl: `Python` · [vaquerizaslab/fanc](https://github.com/vaquerizaslab/fanc) · `GPL-3.0`
  - Existing Rust: none.
  - Priority: `P2`
  - Notes: Niche; useful for differential Hi-C analysis. Not a
    porting target.

- [ ] **`chromosight`** — pattern-matching loop / stripe caller.
  - Reference impl: `Python` · [koszullab/chromosight](https://github.com/koszullab/chromosight) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Algorithm is a 2D convolution / template matching on the
    contact map. `ndarray` + an FFT crate covers it; small focused
    Rust port worth ~100× speedups on dense matrices.

- [ ] **`mustache`** — scale-space chromatin loop caller.
  - Reference impl: `Python` · [ay-lab/mustache](https://github.com/ay-lab/mustache) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Scale-space blob detection on the contact matrix. Like
    chromosight, well suited to `ndarray` + Gaussian-pyramid Rust
    crates.
