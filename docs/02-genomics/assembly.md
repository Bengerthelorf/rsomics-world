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
  - Existing Rust: none for full SPAdes; building-block crates `debruijn`, `ggcat`
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream limited SSE
  - Quadrant: —
  - GPU-amenable: maybe — k-mer counting prequel is GPU-friendly; graph traversal is not
  - Upstream license: `GPL-2.0`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-spades`)
  - Consumes primitives: `rsomics-kmer`, `debruijn` / `ggcat` (foundation), `noodles-fastq`, `noodles-fasta`
  - Notes: GPL licence limits derivative work — clean-room only. A Rust assembler in the SPAdes mold (multi-k de Bruijn + scaffolding) is a Phase-3+ project. `GGCAT` already covers the cDBG core.

- [ ] **`ABySS`** — Bloom-filter-based parallel assembler.
  - Reference impl: `C++` · [bcgsc/abyss](https://github.com/bcgsc/abyss) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream MPI + pthreads
  - SIMD: upstream limited
  - Quadrant: —
  - GPU-amenable: maybe — Bloom-filter lookups are SIMT-trivial
  - Upstream license: `GPL-3.0`
  - Priority: `P2`
  - Layer: `B` (tool — `rsomics-abyss`)
  - Consumes primitives: `fastbloom` ([01-foundations/data-structures](../01-foundations/data-structures.md)), `rsomics-kmer`, `noodles-fastq`, `noodles-fasta`
  - Notes: Distinctive Bloom-filter unitig stage maps well to Rust's `fastbloom` ecosystem. GPL-3.0 forces clean-room work. Lower priority than hifiasm-class work.

- [ ] **`MaSuRCA`** — hybrid short+long-read assembler with super-read technique.
  - Reference impl: `C/C++/Perl` · [alekseyzimin/masurca](https://github.com/alekseyzimin/masurca) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads + Perl orchestration
  - SIMD: upstream limited
  - Quadrant: —
  - GPU-amenable: maybe — only the SW extension steps; pipeline glue is CPU
  - Upstream license: `GPL-3.0`
  - Priority: `P2`
  - Layer: `B` (tool — `rsomics-masurca` if rebuilt)
  - Consumes primitives: `rsomics-kmer`, future `rsomics-chain`, `noodles-bam`
  - Notes: Build system + Perl glue is most of the complexity. Hybrid workflows are increasingly replaced by long-read-only with HiFi.

- [x] **`sparrowhawk`** — bacterial dBG assembler (Rust-native).
  - Reference impl: `Rust` · [bacpop/sparrowhawk](https://github.com/bacpop/sparrowhawk) · `Apache-2.0`
  - Existing Rust: [`sparrowhawk`](https://github.com/bacpop/sparrowhawk) (binary tool, install from source — crates.io name is squatted by an unrelated Shogi library)
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: —
  - Parallelism: rayon
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: maybe — bacterial-scale data fits in GPU memory; engineering cost not justified at this scale
  - Upstream license: `Apache-2.0`
  - Priority: `P2`
  - Layer: `B` (tool — adopt as `rsomics-sparrowhawk`)
  - Consumes primitives: `rsomics-kmer`, `debruijn` or `ggcat`
  - Notes: Niche (bacterial-scale). Treat as a reference for *how* a pure-Rust assembler is structured. Recently updated (2026-05).

### Long-read assemblers

- [ ] **`hifiasm`** — haplotype-resolved HiFi assembler.
  - Reference impl: `C` · [chhylp123/hifiasm](https://github.com/chhylp123/hifiasm) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads, scales well to ~64 threads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — string-graph overlap is irregular; phasing pass is more amenable
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-hifiasm`)
  - Consumes primitives: `rsomics-kmer`, `block-aligner`, `noodles-fasta`, future `rsomics-overlap`
  - Notes: Production standard for diploid HiFi assembly. Pure-Rust port is large but high-value. Phasing-aware overlap graph is the differentiator vs. older OLC assemblers.

- [ ] **`Flye`** — long-read assembler (ONT + PacBio).
  - Reference impl: `C++ / Python` · [mikolmogorov/Flye](https://github.com/mikolmogorov/Flye) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads + Python orchestration
  - SIMD: upstream SSE in alignment kernel
  - Quadrant: —
  - GPU-amenable: maybe — repeat-graph traversal is irregular
  - Upstream license: `BSD-3-Clause`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-flye`)
  - Consumes primitives: `rsomics-kmer`, `block-aligner`, `noodles-fasta`, future `rsomics-repeat-graph`
  - Notes: ONT default since the R9 era. Uses repeat-graph approach different from hifiasm's string-graph. BSD licence is friendly.

- [ ] **`Verkko`** — T2T-class hybrid HiFi + ONT assembler.
  - Reference impl: `C++ / Python / Snakemake` · [marbl/verkko](https://github.com/marbl/verkko) · `Public domain / BSD`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Snakemake pipeline orchestrates hifiasm + MBG + GraphAligner
  - SIMD: inherited from constituent tools
  - Quadrant: —
  - GPU-amenable: maybe — each constituent has its own profile
  - Upstream license: `Public domain` (core); `BSD-3-Clause` parts
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-verkko` as a Snakemake-equivalent wrapper)
  - Consumes primitives: future `rsomics-hifiasm`, `rsomics-graph-align`, workflow engine ([09-workflow-utility](../09-workflow-utility/workflow-engines.md))
  - Notes: The T2T human reference assembly tool. Heavy pipeline. Realistic Rust work is wrapping the DAG, not rewriting the core.

- [ ] **`Canu`** — OLC long-read assembler (legacy).
  - Reference impl: `C++/Perl` · [marbl/canu](https://github.com/marbl/canu) · `GPL-3.0` (dual)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads + grid engines
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — OLC is irregular; not a fit
  - Upstream license: `GPL-3.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Authors declare end-of-life (2021) and recommend Flye/hifiasm/Verkko instead. Document only.

- [ ] **`wtdbg2`** — fuzzy-Bruijn-graph long-read assembler.
  - Reference impl: `C` · [ruanjue/wtdbg2](https://github.com/ruanjue/wtdbg2) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — fuzzy-Bruijn graph; same constraints as other dBG variants
  - Upstream license: `GPL-3.0`
  - Priority: `P2`
  - Layer: `B` (tool — `rsomics-wtdbg2`)
  - Consumes primitives: `rsomics-kmer`, `block-aligner`
  - Notes: Fast but lower contiguity than hifiasm/Flye. Used in some niche large-genome pipelines.

- [ ] **`NextDenovo`** — long-read assembler from Nextomics.
  - Reference impl: `C/Python` · [Nextomics/NextDenovo](https://github.com/Nextomics/NextDenovo) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads + Python orchestration
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe
  - Upstream license: `GPL-3.0`
  - Priority: `P2`
  - Layer: `B` (tool — `rsomics-nextdenovo`)
  - Consumes primitives: `rsomics-kmer`, `block-aligner`
  - Notes: Used in plant + animal genome projects, especially in China. GPL licence.

- [ ] **`Shasta`** — fast Nanopore assembler.
  - Reference impl: `C++` · [paoloshasta/shasta](https://github.com/paoloshasta/shasta) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads (memory-mapped index)
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — memory-bound index probing
  - Upstream license: `BSD-3-Clause`
  - Priority: `P2`
  - Layer: `B` (tool — `rsomics-shasta`)
  - Consumes primitives: `rsomics-kmer`, `block-aligner`, `noodles-fasta`
  - Notes: Targets large-scale clinical ONT. Memory-hungry. Less common in 2026 pipelines than Flye.

- [ ] **`Raven`** — fast long-read assembler.
  - Reference impl: `C++` · [lbcb-sci/raven](https://github.com/lbcb-sci/raven) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe
  - Upstream license: `MIT`
  - Priority: `P2`
  - Layer: `B` (tool — `rsomics-raven`)
  - Consumes primitives: `rsomics-kmer`, `block-aligner`, `noodles-fasta`
  - Notes: Small codebase, MIT licence — a tractable Rust port if a contributor wants a benchmark target.

- [x] **`rust-mdbg`** — minimizer-space dBG long-read assembler (Rust-native).
  - Reference impl: `Rust` · [ekimb/rust-mdbg](https://github.com/ekimb/rust-mdbg) · `MIT`
  - Existing Rust: [`rust-mdbg`](https://github.com/ekimb/rust-mdbg) (binary tool, install from source); companion library [`rust-seq2kminmers`](https://crates.io/crates/rust-seq2kminmers) `0.1.0`
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: —
  - Parallelism: rayon
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: maybe — minimizer extraction is SIMT-friendly; the assembly stage is not
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `B` (tool — adopt as `rsomics-mdbg`)
  - Consumes primitives: `rsomics-kmer`, `rust-seq2kminmers`
  - Notes: 52× human genome HiFi assembly in ~10 min on 8 threads. Research-grade quality (less polished than hifiasm) but proves a pure-Rust HiFi assembler is feasible. Adopt and extend.

### Polishers

- [ ] **`Pilon`** — short-read-based polisher.
  - Reference impl: `Scala/Java` · [broadinstitute/pilon](https://github.com/broadinstitute/pilon) · `GPL-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: JVM threading
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — primarily pileup walking
  - Upstream license: `GPL-2.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Slow, memory-hungry, largely superseded by NextPolish / Polypolish for modern workflows. Document only.

- [ ] **`Racon`** — long-read consensus polisher.
  - Reference impl: `C++` · [lbcb-sci/racon](https://github.com/lbcb-sci/racon) · `MIT`
  - Existing Rust: [`rust-spoa`](https://crates.io/crates/rust-spoa) `0.2.4` covers the partial-order alignment building block (FFI wrapper of the C++ SPOA library)
  - Existing Rust kind: `FFI-wrapper` (rust-spoa); `none` for full Racon
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: ② (rust-spoa); — (full Racon)
  - GPU-amenable: maybe — SPOA has GPU variants in the literature
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-racon`)
  - Consumes primitives: future `rsomics-spoa` (pure-Rust SPOA), `noodles-fasta`, `block-aligner`
  - Notes: Heavy partial-order alignment. `rust-spoa` is FFI to the same SPOA library Racon uses — a pure-Rust SPOA + Racon port is achievable.

- [ ] **`Medaka`** — neural-net ONT polisher.
  - Reference impl: `Python (TF/Keras)` · [nanoporetech/medaka](https://github.com/nanoporetech/medaka) · `MPL-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: TF/Keras GPU
  - SIMD: TF/Keras kernels
  - Quadrant: —
  - GPU-amenable: yes — neural-net inference, dense ops map directly to SIMT
  - Upstream license: `MPL-2.0`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-medaka`)
  - Consumes primitives: `candle` or `burn` for inference, `noodles-fasta` / `noodles-bam`, future `rsomics-pileup`
  - Notes: Inference workload — port the model to `candle` / `burn` rather than re-train. Tracks closely with PEPPER-Margin-DeepVariant in [`variant-calling.md`](variant-calling.md).

- [ ] **`NextPolish`** — short-read polisher.
  - Reference impl: `C/Python` · [Nextomics/NextPolish](https://github.com/Nextomics/NextPolish) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `Polypolish` (C++)
  - Parallelism: upstream pthreads + Python orchestration
  - SIMD: upstream limited
  - Quadrant: —
  - GPU-amenable: no — alignment-based pileup walking
  - Upstream license: `BSD-3-Clause`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-polish` covering both short-read polishing modes)
  - Consumes primitives: `noodles-bam`, future `rsomics-pileup`
  - Notes: Reportedly highest accuracy among short-read polishers; pairs with Medaka in the most-common modern polishing combo. BSD licence is friendly.
