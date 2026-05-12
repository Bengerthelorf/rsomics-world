# Assembly

> *De novo* genome assembly: from short and long reads to contigs / scaffolds,
> plus polishing.

## Scope

Short-read assemblers (SPAdes, ABySS, MaSuRCA), long-read assemblers
(Flye, hifiasm, Canu, wtdbg2, NextDenovo, Shasta, Raven, Verkko), and
polishers (Pilon, Racon, Medaka, NextPolish). Excludes metagenomic
assembly (module 06), transcriptome assembly (module 03), and pangenome
graph construction (module 08). De Bruijn graph *data structures* live in
[`01-foundations/data-structures.md`](../01-foundations/data-structures.md);
this doc is about the full assembler pipelines built on top.

## Design notes

- The assembler space has moved dramatically toward long reads
  (HiFi/ONT). For HiFi, **hifiasm** is the production standard;
  **Verkko** for telomere-to-telomere assemblies combining HiFi + ONT.
  Canu has been declared end-of-life by its authors (2021).
- Pure-Rust assembly *is* possible — see `rust-mdbg` (minimizer-space
  dBG for HiFi) and `GGCAT` (compacted/coloured dBG) — but no Rust
  assembler approaches the feature completeness of hifiasm or Flye yet.
- Short-read assembly is largely a solved-and-stagnant problem; SPAdes
  remains the default but rarely sees major performance work.
- Polishers are increasingly important as long-read accuracy improves:
  Medaka (ONT, neural-net) and NextPolish (Illumina, alignment-based)
  are the modern pair. Pilon is legacy but still cited.
- Memory-mapping and BGZF-friendly IO matter at assembly scale (10s of
  TB intermediate data for large vertebrates). The foundations layer
  (`noodles-bgzf`, `niffler`) carries most of the weight here.

## TODO

### Short-read assemblers

- [ ] **`SPAdes`** — multi-mode short-read / hybrid assembler.
  - Reference impl: `C++` · [ablab/spades](https://github.com/ablab/spades) · `GPL-2.0`
  - Existing Rust: none verified for full SPAdes; building-block crates
    `rust-debruijn`, `ggcat`
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: GPL licence limits derivative work. A Rust assembler in
    the SPAdes mold (multi-k de Bruijn + scaffolding) is a Phase-3+
    project. `GGCAT` already covers the cDBG core.

- [ ] **`ABySS`** — Bloom-filter-based parallel assembler.
  - Reference impl: `C++` · [bcgsc/abyss](https://github.com/bcgsc/abyss) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Distinctive Bloom-filter unitig stage maps well to Rust's
    `fastbloom` ecosystem. GPL-3.0 again forces clean-room work. Lower
    priority than hifiasm-class work.

- [ ] **`MaSuRCA`** — hybrid short+long-read assembler with super-read
  technique.
  - Reference impl: `C/C++/Perl` · [alekseyzimin/masurca](https://github.com/alekseyzimin/masurca) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Build system + Perl glue is most of the complexity. Hybrid
    workflows are increasingly replaced by long-read-only with HiFi.

- [ ] **`sparrowhawk`** — bacterial dBG assembler (Rust, new).
  - Reference impl: `Rust` · [bacpop/sparrowhawk](https://github.com/bacpop/sparrowhawk) · `Apache-2.0`
  - Existing Rust: [`sparrowhawk`](https://github.com/bacpop/sparrowhawk) (itself)
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Niche (bacterial-scale). Treat as a reference for *how* a
    pure-Rust assembler is structured.

### Long-read assemblers

- [ ] **`hifiasm`** — haplotype-resolved HiFi assembler.
  - Reference impl: `C` · [chhylp123/hifiasm](https://github.com/chhylp123/hifiasm) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Production standard for diploid HiFi assembly. Pure-Rust port
    is large but high-value. Phasing-aware overlap graph is the
    differentiator vs. older OLC assemblers.

- [ ] **`Flye`** — long-read assembler (ONT + PacBio).
  - Reference impl: `C++ / Python` · [mikolmogorov/Flye](https://github.com/mikolmogorov/Flye) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: ONT default since the R9 era. Uses repeat-graph approach
    different from hifiasm's string-graph. BSD licence is friendly.

- [ ] **`Verkko`** — T2T-class hybrid HiFi + ONT assembler.
  - Reference impl: `C++ / Python / Snakemake` · [marbl/verkko](https://github.com/marbl/verkko) · `Public domain / BSD`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: The T2T human reference assembly tool. Heavy pipeline (calls
    into hifiasm, MBG, GraphAligner). Realistic Rust work is wrapping
    the Snakemake DAG, not rewriting the core.

- [ ] **`Canu`** — OLC long-read assembler (legacy).
  - Reference impl: `C++/Perl` · [marbl/canu](https://github.com/marbl/canu) · `GPL-3.0` (dual)
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Authors declare end-of-life (2021) and recommend Flye/hifiasm/
    Verkko instead. Document only.

- [ ] **`wtdbg2`** — fuzzy-Bruijn-graph long-read assembler.
  - Reference impl: `C` · [ruanjue/wtdbg2](https://github.com/ruanjue/wtdbg2) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Fast but lower contiguity than hifiasm/Flye. Used in some
    niche large-genome pipelines.

- [ ] **`NextDenovo`** — long-read assembler from Nextomics.
  - Reference impl: `C/Python` · [Nextomics/NextDenovo](https://github.com/Nextomics/NextDenovo) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Used in plant + animal genome projects, especially in China.
    GPL licence.

- [ ] **`Shasta`** — fast Nanopore assembler.
  - Reference impl: `C++` · [paoloshasta/shasta](https://github.com/paoloshasta/shasta) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Targets large-scale clinical ONT. Memory-hungry. Less common
    in 2026 pipelines than Flye.

- [ ] **`Raven`** — fast long-read assembler.
  - Reference impl: `C++` · [lbcb-sci/raven](https://github.com/lbcb-sci/raven) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Small codebase, MIT licence — a tractable Rust port if a
    contributor wants a benchmark target.

- [~] **`rust-mdbg`** — minimizer-space dBG long-read assembler (Rust).
  - Reference impl: `Rust` · [ekimb/rust-mdbg](https://github.com/ekimb/rust-mdbg) · `MIT`
  - Existing Rust: [`rust-mdbg`](https://github.com/ekimb/rust-mdbg)
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: 52× human genome HiFi assembly in ~10 min on 8 threads.
    Research-grade quality (less polished than hifiasm) but proves a
    pure-Rust HiFi assembler is feasible. Adopt and extend.

### Polishers

- [ ] **`Pilon`** — short-read-based polisher.
  - Reference impl: `Scala/Java` · [broadinstitute/pilon](https://github.com/broadinstitute/pilon) · `GPL-2.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Slow, memory-hungry, largely superseded by NextPolish /
    Polypolish for modern workflows. Document only.

- [ ] **`Racon`** — long-read consensus polisher.
  - Reference impl: `C++` · [lbcb-sci/racon](https://github.com/lbcb-sci/racon) · `MIT`
  - Existing Rust: [`rust-spoa`](https://crates.io/crates/rust-spoa)
    covers the partial-order alignment building block
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Heavy partial-order alignment. `rust-spoa` is FFI to the same
    SPOA library Racon uses — a pure-Rust SPOA + Racon port is
    achievable.

- [ ] **`Medaka`** — neural-net ONT polisher.
  - Reference impl: `Python (TF/Keras)` · [nanoporetech/medaka](https://github.com/nanoporetech/medaka) · `MPL-2.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Inference workload — port the model to `candle` / `burn` rather
    than re-train. Tracks closely with PEPPER-Margin-DeepVariant in
    [`variant-calling.md`](variant-calling.md).

- [ ] **`NextPolish`** — short-read polisher.
  - Reference impl: `C/Python` · [Nextomics/NextPolish](https://github.com/Nextomics/NextPolish) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing non-C alternatives: `Polypolish` (C++)
  - Priority: `P1`
  - Notes: Reportedly highest accuracy among short-read polishers; pairs
    with Medaka in the most-common modern polishing combo. BSD licence is
    friendly.
