# Multiple sequence alignment

> Progressive, iterative, consistency-based, and HMM-divide-and-conquer
> multiple sequence aligners for DNA, RNA, and protein.

## Scope

Includes: classical pairwise-progressive (Clustal Omega, T-Coffee),
iterative refinement (MUSCLE5, MAFFT), ultra-scale linear-time aligners
(FAMSA, KAlign3), and divide-and-conquer for highly diverged / ultra-large
inputs (UPP, PASTA). Excludes: pairwise alignment (covered under
foundations / genomics) and structure-based MSA (Foldseek family, see
[../07-proteomics-structure/structure-analysis.md](../07-proteomics-structure/structure-analysis.md)).

## Design notes

- MAFFT, MUSCLE5, and Clustal Omega are mature, well-optimized C/C++
  with stable APIs. Marginal benefit of a Rust port is mainly licensing
  and embeddability — none of them are LGPL/MIT, so calling them as
  libraries from a permissive Rust application is awkward.
- FAMSA is the fastest at ultra-scale (millions of sequences) and ships
  AVX2/AVX-512 inner loops. A pure-Rust FAMSA is interesting but the
  SIMD work is the entire point.
- KAlign3 is small (~6k LOC), fast, and BSD-2 — the cleanest pure-Rust
  port target in this list.
- T-Coffee's value is consistency-based scoring and meta-alignment
  combination; a Rust analogue could be small and elegant.
- UPP and PASTA are Python orchestrators that wrap HMMER + an MSA tool +
  a tree tool. The Python is replaceable; the underlying algorithms
  (SEPP-style decomposition) are interesting.
- License watch: MAFFT **BSD-3**, MUSCLE5 **GPL-3**, Clustal Omega
  **GPL-2**, T-Coffee **GPL-3**, KAlign3 **BSD-2-Clause**, FAMSA **GPL-3**,
  UPP/PASTA mixed (Python GPL-3-ish).

## TODO

- [ ] **`MAFFT`** — fast, accurate iterative MSA with multiple modes (FFT-NS, L-INS-i, etc.).
  - Reference impl: `C` · [GSLBiotech/mafft](https://github.com/GSLBiotech/mafft) · `BSD-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE/AVX
  - Quadrant: —
  - GPU-amenable: maybe — DP-style alignment scoring SIMT-friendly
  - Upstream license: `BSD-3`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-mafft`)
  - Consumes primitives: `noodles-fasta`, `block-aligner`, `nalgebra`, future `rsomics-hmm`
  - Notes: De-facto MSA workhorse. Large codebase (~50k LOC). Rust port is a multi-month project. Lower priority than FAMSA/KAlign3 since MAFFT is BSD-3 and already integrates cleanly via FFI for the near term.

- [ ] **`MUSCLE5`** — ensemble HMM alignment with permuted guide trees.
  - Reference impl: `C++` · [rcedgar/muscle](https://github.com/rcedgar/muscle) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE/AVX
  - Quadrant: —
  - GPU-amenable: maybe — ensemble alignment parallelises trivially
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-msa` (ensemble mode)
  - Consumes primitives: `noodles-fasta`, `block-aligner`, future `rsomics-hmm`
  - Notes: Algorithmically distinct (ensemble alignment). GPL-3 makes direct embedding into permissive Rust apps awkward. Pure-Rust port would need to reproduce the published ensemble scheme — non-trivial.

- [ ] **`Clustal Omega`** — HMM-based MSA, fast for large alignments.
  - Reference impl: `C` · [clustal.org/omega](http://www.clustal.org/omega/) · `GPL-2`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — HMM scoring
  - Upstream license: `GPL-2`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Older codebase, slow community development. Skip; MAFFT/FAMSA cover this niche better.

- [ ] **`T-Coffee`** — consistency-based MSA and meta-alignment.
  - Reference impl: `C` · [cbcrg/tcoffee](https://github.com/cbcrg/tcoffee) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — consistency scoring per pair
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-msa` (consistency / meta-alignment mode)
  - Consumes primitives: `block-aligner`, `noodles-fasta`
  - Notes: Niche but algorithmically interesting (consistency scoring). Small scoped Rust port could be valuable for embedding in `rsomics-msa` as a meta-aligner.

- [ ] **`KAlign3`** — fast linear-time aligner via LCS heuristic.
  - Reference impl: `C` · [TimoLassmann/kalign](https://github.com/TimoLassmann/kalign) · `BSD-2-Clause`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — LCS heuristic per-pair
  - Upstream license: `BSD-2-Clause`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-kalign`)
  - Consumes primitives: `noodles-fasta`, `block-aligner`
  - Notes: Best pure-Rust port target in this list: small (~6k LOC), permissive license, fast, modern algorithm. `rsomics-kalign` ships as a small crate with a clean CLI. Excellent starter project.

- [ ] **`FAMSA`** — ultra-scale progressive MSA with SIMD-tuned kernels.
  - Reference impl: `C++` · [refresh-bio/FAMSA](https://github.com/refresh-bio/FAMSA) · `GPL-3`
  - Existing Rust: none verified. `pyfamsa` ([althonos/pyfamsa](https://github.com/althonos/pyfamsa)) provides Python+Cython bindings only
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream AVX2/AVX-512 hand intrinsics
  - Quadrant: —
  - GPU-amenable: yes — DP kernels SIMT-friendly
  - Upstream license: `GPL-3`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-famsa`)
  - Consumes primitives: `noodles-fasta`, `block-aligner`, `wide`/`pulp` for SIMD, `rayon`
  - Notes: Best for million-sequence alignments (e.g. Pfam-scale). The SIMD inner loops are the value; Rust needs `std::simd` or `wide`/`pulp` to be competitive. GPL-3 = clean-room rewrite.

- [ ] **`UPP`** — SEPP-based MSA for ultra-large / fragmentary datasets.
  - Reference impl: `Python` · [smirarab/sepp](https://github.com/smirarab/sepp) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: inherited
  - Quadrant: —
  - GPU-amenable: maybe — divide-and-conquer over fragments parallelises
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-msa` (UPP / SEPP mode)
  - Consumes primitives: future `rsomics-hmm`, `rsomics-kalign` or `rsomics-famsa`, `noodles-fasta`
  - Notes: Python orchestrator; meaningful logic in HMMER + an MSA tool underneath. Port only after the underlying tools exist.

- [ ] **`PASTA`** — practical SATe-aligner for large datasets.
  - Reference impl: `Python` · [smirarab/PASTA](https://github.com/smirarab/PASTA) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: inherited
  - Quadrant: —
  - GPU-amenable: maybe — same family as UPP
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-msa` (PASTA / SATe mode)
  - Consumes primitives: same as UPP
  - Notes: Same shape as UPP. Defer.
