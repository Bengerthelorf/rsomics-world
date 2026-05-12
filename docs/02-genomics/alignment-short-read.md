# Short-read alignment

> Aligning Illumina-class reads (50–300 bp, low error rate) to a reference
> genome.

## Scope

Aligners specialised for short, accurate reads — BWA family, Bowtie2,
SNAP, Strobealign, NovoAlign. Excludes long-read aligners (see
[`alignment-long-read.md`](alignment-long-read.md)), pseudo-alignment for
quantification (module 03), and graph genome alignment (treated as a
separate problem under module 08).

## Design notes

- BWA-MEM is still the *de-facto* short-read aligner for clinical
  pipelines in 2026. **bwa-mem2** (the modern C++/SIMD rewrite from the
  same group) is 2–3× faster on AVX-512 hardware and is the right
  performance baseline for any Rust rewrite to compare against — not
  the original 2013-era `bwa`.
- The inner loop is banded Smith-Waterman extension. SIMD (SSE4 / AVX2 /
  AVX-512 / NEON) is the largest single performance lever; pure-Rust
  ports must use `std::simd` or `pulp` to remain competitive.
- The FM-index for the reference is rebuildable in minutes for human
  but takes >1 hour for plant genomes — a parallel suffix-array
  constructor (see
  [`01-foundations/data-structures.md`](../01-foundations/data-structures.md))
  is on the critical path.
- Strobealign is the most credible *new* short-read aligner since
  bwa-mem2: dynamic seed sizes via syncmer-thinned strobemers. C++,
  permissive licence, already compiles against Rust components.
- NovoAlign remains commercial / closed-source. We document it for
  completeness but cannot port or wrap it.
- 10x Genomics ships [`rust-bwa`](https://github.com/10XGenomics/rust-bwa) —
  FFI only; no pure-Rust BWA implementation exists.

## TODO

- [~] **`bwa-mem` / `bwa-mem2`** — seed-and-extend Burrows-Wheeler aligner.
  - Reference impl: `C` · [lh3/bwa](https://github.com/lh3/bwa) · `MIT/GPL-3.0` (dual)
  - Modern reference impl: `C++` · [bwa-mem2/bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) · `MIT`
  - Existing Rust: [`rust-bwa`](https://github.com/10XGenomics/rust-bwa) (FFI wrapper, BWA C API);
    no pure-Rust port verified
  - Existing non-C alternatives: `bwa-meme` (learned-index C++ fork);
    `BWA-MEME` uses Rust for *training* the learned index, not alignment
  - Priority: `P0`
  - Notes: Inner SW kernel is SIMD-critical. Start with the FFI wrapper
    to unblock downstream pipelines; plan pure-Rust port after the
    foundations FM-index work lands. Compare against `bwa-mem2` not
    `bwa 0.7.17` for fairness.

- [ ] **`Bowtie2`** — gapped seed-extend aligner with end-to-end and local
  modes.
  - Reference impl: `C++` · [BenLangmead/bowtie2](https://github.com/BenLangmead/bowtie2) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: GPL-3.0 license complicates re-derivation; a clean-room Rust
    port can ship under MIT/Apache-2.0 as a sibling. Bowtie2 retains a
    loyal user base (epigenetics, ATAC-seq) so module 05 will need it.

- [ ] **`SNAP`** — hash-based aligner, parallel-friendly.
  - Reference impl: `C++` · [amplab/snap](https://github.com/amplab/snap) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Apache-2.0 licence is friendly. Hash-table index simplifies
    parallelism vs. FM-index. Strong scaling claims but limited adoption
    in 2026 pipelines — keep on the radar but not a priority.

- [ ] **`Strobealign`** — strobemer seed-extend aligner.
  - Reference impl: `C++` · [ksahlin/strobealign](https://github.com/ksahlin/strobealign) · `MIT`
  - Existing Rust: none verified (strobealign itself requires Rust at *build*
    time for an auxiliary component, but the aligner is C++)
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Faster than bwa-mem2 on 150–300 bp reads with comparable
    accuracy. Algorithm is younger (less battle-tested) but the
    code is small and re-derivable. Strong candidate for an *early*
    Rust rewrite where we are not chasing a moving target.

- [ ] **`NovoAlign`** — commercial proprietary aligner.
  - Reference impl: closed-source · [novocraft.com](https://www.novocraft.com/) · proprietary
  - Existing Rust: n/a
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Cannot port (closed source). Document interop only:
    `noodles-sam` should parse its outputs without warnings.

- [~] **`block-aligner`** — SIMD-accelerated banded SW kernel (reusable
  by any new aligner).
  - Reference impl: `Rust` · [Daniel-Liu-c0deb0t/block-aligner](https://github.com/Daniel-Liu-c0deb0t/block-aligner) · `MIT`
  - Existing Rust: [`block-aligner`](https://crates.io/crates/block-aligner)
  - Existing non-C alternatives: `ksw2` (C, used by minimap2); `WFA2-lib` (C)
  - Priority: `P0`
  - Notes: Adopt as the standard SW kernel for any pure-Rust short-read
    aligner work. SSE2/AVX2/NEON/WASM-SIMD paths already shipped. Avoids
    re-implementing the most-vectorised inner loop in the field.
