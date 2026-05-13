# Structural-variant calling

> Large (≥50 bp) variant callers — deletions, duplications, inversions,
> translocations, insertions, CNVs.

## Scope

Short-read SV callers (Manta, Delly, Lumpy, SvABA, GRIDSS, Wham) and
long-read SV callers (Sniffles, cuteSV, SVIM, pbsv, Sawfish). Small
variants (<50 bp) live in [`variant-calling.md`](variant-calling.md);
copy-number analysis on arrays / WGS coverage is a sibling concern that
will live in module 05 or here, TBD.

## Design notes

- Short-read SV callers rely on three signals: read pairs (insert-size
  anomalies), split reads (clipped alignments), and local assembly.
  GRIDSS/Manta/SvABA combine all three; Lumpy/Delly skip assembly.
  Consensus calling (typically using `survivor` or `Jasmine`) often
  outperforms any single caller.
- Long-read SV callers exploit alignment-spanning reads — much simpler
  algorithmically than short-read. Sniffles + cuteSV are the production
  standards in 2026; both have minimap2-based input pipelines.
- Sawfish (PacBio) and Severus (Kolmogorov lab) are the newest entrants
  worth tracking — joint SV+CNV calling for HiFi, somatic SV calling for
  long reads.
- A pure-Rust short-read SV caller is a large undertaking but tractable
  in stages: pair-end + split-read signal extraction is straightforward
  with `noodles-bam`; local assembly leverages
  [`01-foundations/data-structures.md`](../01-foundations/data-structures.md)
  cDBG work.
- For long-read SV calling, the bottleneck is well-defined: read
  parsing + cluster + filter. A `rsomics-sniffles` could land before a
  `rsomics-manta`.

## TODO

### Short-read SV callers

- [ ] **`Manta`** — pair-end + split-read + local-assembly SV caller.
  - Reference impl: `C++` · [Illumina/manta](https://github.com/Illumina/manta) · `GPL-3.0` — **upstream repo archived**
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE in local assembly
  - Quadrant: —
  - GPU-amenable: maybe — pair/split-read filtering is SIMT-friendly; assembly less so
  - Upstream license: `GPL-3.0`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-manta`)
  - Consumes primitives: `noodles-bam`, `noodles-vcf`, `rsomics-kmer`, `debruijn` / `ggcat`, `block-aligner`
  - Notes: Most-used short-read SV caller. Upstream is archived — Illumina is no longer maintaining; ecosystem still uses it heavily. GPL-3.0 licence; clean-room re-implementation required. Strong DEL/DUP/INS/BND coverage, moderate INV. Used by GATK-SV pipeline as one of three callers.

- [ ] **`Delly`** — pair-end + split-read SV caller.
  - Reference impl: `C++` · [dellytools/delly](https://github.com/dellytools/delly) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — signal extraction is parallel; clustering less so
  - Upstream license: `BSD-3-Clause`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-delly`)
  - Consumes primitives: `noodles-bam`, `noodles-vcf`, `rsomics-intervals`
  - Notes: BSD-3-Clause is friendly. Smaller codebase than Manta. Strong on DEL/DUP/INV/BND. Often run alongside Manta for consensus.

- [ ] **`Lumpy`** — probabilistic SV discovery (legacy).
  - Reference impl: `C++` · [arq5x/lumpy-sv](https://github.com/arq5x/lumpy-sv) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `smoove` (Go wrapper for the same pipeline)
  - Parallelism: upstream pthreads via samtools/bcftools workers
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — probabilistic clustering, irregular
  - Upstream license: `MIT`
  - Priority: `P2`
  - Layer: `B` (tool — `rsomics-lumpy`)
  - Consumes primitives: `noodles-bam`, `noodles-vcf`, `rsomics-stats`
  - Notes: Largely superseded by Manta + Delly. MIT licence and clean architecture make it a tractable benchmark target.

- [ ] **`SvABA`** — local-assembly somatic + germline SV caller.
  - Reference impl: `C++` · [walaj/svaba](https://github.com/walaj/svaba) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE in SW
  - Quadrant: —
  - GPU-amenable: maybe — SW kernel is SIMT-amenable
  - Upstream license: `GPL-3.0`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-svaba`)
  - Consumes primitives: `noodles-bam`, `noodles-vcf`, `debruijn`, `block-aligner`
  - Notes: Best-in-class for small SVs and short insertions in cancer panels. Heavy local-assembly path. Pairs naturally with our `debruijn` / `block-aligner` work.

- [ ] **`GRIDSS`** — graph-based break-end caller.
  - Reference impl: `Java` · [PapenfussLab/gridss](https://github.com/PapenfussLab/gridss) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: JVM
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — graph traversal irregular, but signal extraction is parallel
  - Upstream license: `GPL-3.0`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-gridss`)
  - Consumes primitives: `noodles-bam`, `noodles-vcf`, `debruijn`, `block-aligner`
  - Notes: Highest sensitivity among short-read callers but slowest and most resource-hungry. Outputs raw BND records; downstream tools (LINX, GRIDSS-PURPLE-LINX) categorise into typed events.

- [ ] **`Wham`** — split-read + pair SV caller.
  - Reference impl: `C++` · [zeeev/wham](https://github.com/zeeev/wham) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — clustering-dominated
  - Upstream license: `MIT`
  - Priority: `P2`
  - Layer: `B` (tool — `rsomics-wham`)
  - Consumes primitives: `noodles-bam`, `noodles-vcf`
  - Notes: Used in GATK-SV pipeline. MIT licence. Maintenance is slow upstream. Lower priority than Manta/Delly.

### Long-read SV callers

- [ ] **`Sniffles`** — long-read SV caller (PacBio + ONT).
  - Reference impl: `Python / C++` · [fritzsedlazeck/Sniffles](https://github.com/fritzsedlazeck/Sniffles) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream Python threading + pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: no — clustering / filtering, latency-dominated
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-sniffles`)
  - Consumes primitives: `noodles-bam`, `noodles-vcf`, `rsomics-intervals`
  - Notes: Production standard. Sniffles2 (current) is Python-heavy and relatively small (~5k LOC core). MIT licence. High-priority Rust target — produces VCFs that the rest of the pipeline (Jasmine, Truvari) accepts natively.

- [ ] **`cuteSV`** — long-read SV caller with high precision/recall.
  - Reference impl: `Python` · [tjiangHIT/cuteSV](https://github.com/tjiangHIT/cuteSV) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-cutesv`)
  - Consumes primitives: `noodles-bam`, `noodles-vcf`, `rsomics-intervals`
  - Notes: Often ranks at or near top in benchmarks (sometimes ahead of Sniffles). MIT, pure Python — a Rust port should be straightforward and offer large speedups (Python is the bottleneck, not the algorithm).

- [ ] **`SVIM`** — long-read SV caller (Python).
  - Reference impl: `Python` · [eldariont/svim](https://github.com/eldariont/svim) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python single-threaded
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no
  - Upstream license: `GPL-3.0`
  - Priority: `P2`
  - Layer: `B` (tool — `rsomics-svim`)
  - Consumes primitives: `noodles-bam`, `noodles-vcf`
  - Notes: Slower than Sniffles/cuteSV; GPL-3.0 complicates derivation. Active research lab maintenance.

- [ ] **`pbsv`** — PacBio's official SV caller.
  - Reference impl: `C++` (closed-source binary release) · [PacificBiosciences/pbsv](https://github.com/PacificBiosciences/pbsv) · `BSD-3-Clause-Clear` (binary-only)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: pthreads (binary)
  - SIMD: not documented
  - Quadrant: —
  - GPU-amenable: unknown — binary-only
  - Upstream license: `BSD-3-Clause-Clear` (binary-only redistribution)
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Binary-only distribution from PacBio. Cannot port directly. PacBio has shifted recommendation to **Sawfish** for HiFi SV/CNV.

- [x] **`Sawfish`** — joint SV + CNV caller for HiFi (PacBio).
  - Reference impl: `Rust` · [PacificBiosciences/sawfish](https://github.com/PacificBiosciences/sawfish) · `BSD-3-Clause-Clear`
  - Existing Rust: [`sawfish`](https://github.com/PacificBiosciences/sawfish) (binary tool, install from source — crates.io name is squatted by `sawfish-client`, an unrelated window-manager crate)
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: —
  - Parallelism: rayon
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: no — clustering + statistical inference, latency-dominated
  - Upstream license: `BSD-3-Clause-Clear`
  - Priority: `P1`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Already Rust. Adopt as the canonical HiFi joint SV/CNV caller. Document interop with our `rsomics-vcf` stack.

- [ ] **`Severus`** — somatic long-read SV caller.
  - Reference impl: `Python` · [KolmogorovLab/Severus](https://github.com/KolmogorovLab/Severus) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no
  - Upstream license: `BSD-3-Clause`
  - Priority: `P2`
  - Layer: `B` (tool — `rsomics-severus`)
  - Consumes primitives: `noodles-bam`, `noodles-vcf`, `rsomics-intervals`
  - Notes: Newer (2024). Same lab as Flye. Worth tracking; Python source is small and tractable.
