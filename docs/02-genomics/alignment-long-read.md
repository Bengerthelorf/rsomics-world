# Long-read alignment

> Aligning PacBio HiFi/CLR and Oxford Nanopore reads (1–100+ kb, error rate
> 0.1%–10%) to a reference genome.

## Scope

Long-read DNA aligners. minimap2 dominates, but several specialised
aligners exist for repetitive references (Winnowmap), SV-aware mapping
(NGMLR), or fast spliced/genomic switching (LRA). RNA-specific long-read
alignment (e.g. IsoQuant, FLAIR) lives in module 03; short-read alignment
in [`alignment-short-read.md`](alignment-short-read.md).

## Design notes

- minimap2 is the gold standard. Its inner alignment kernel (`ksw2`) is
  hand-vectorised SSE/AVX; its chaining algorithm has been rewritten
  several times for performance. A pure-Rust port is an enormous
  undertaking with limited upside *unless* we can integrate it more
  naturally with downstream Rust tools.
- The pragmatic path for the near term: `minimap2-rs` FFI wrapper for
  pipeline integration, with a longer-term goal of replacing the kernel
  with `block-aligner` and the chainer with a pure-Rust port.
- Winnowmap is the right choice for repetitive references (centromeres,
  T2T assembly); NGMLR for SV-heavy nanopore data; LRA for fast all-
  purpose; lordFAST for PacBio CLR (older tech).
- ONT R10.4.1 + duplex base-callers have shifted long-read error
  profiles closer to HiFi over the past two years — many of the
  "ONT-specialised" tools (NGMLR, lordFAST) are losing relevance to
  minimap2 with sensible parameters.

## TODO

- [~] **`minimap2`** — versatile pairwise aligner, the *de-facto* long-read aligner.
  - Reference impl: `C` · [lh3/minimap2](https://github.com/lh3/minimap2) · `MIT`
  - Existing Rust: [`minimap2`](https://crates.io/crates/minimap2) `0.1.31+minimap2.2.30` ([jguhlin/minimap2-rs](https://github.com/jguhlin/minimap2-rs)); paired [`minimap2-sys`](https://crates.io/crates/minimap2-sys) `0.1.30+minimap2.2.30` (raw FFI); [`minimap2-temp`](https://crates.io/crates/minimap2-temp) `0.1.33+minimap2.2.28` (transient pin variant)
  - Existing Rust kind: `FFI-wrapper`
  - Existing non-C alternatives: `mm2-fast` (C++/SIMD fork by bwa-mem2 group)
  - Parallelism: inherits upstream pthreads (`-t` knob)
  - SIMD: inherits upstream's hand-written SSE / AVX intrinsics in ksw2
  - Quadrant: ②
  - GPU-amenable: maybe — NVIDIA Parabricks ships a GPU minimap2; the chainer is the harder half to port
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-minimap` after the kernel and chainer port land; for now adopt FFI)
  - Consumes primitives: `block-aligner` (eventually replaces ksw2), `noodles-bam` for output, future `rsomics-chain`
  - Notes: FFI-only Rust today. Pure-Rust port is the long-term goal but is a multi-person-year effort: chainer, ksw2-like SW, presets, splice model. Adopt `minimap2` crate (`jguhlin/minimap2-rs`) for now; track `mm2-fast` upstream improvements for the SIMD baseline.

- [ ] **`NGMLR`** — long-read mapper specialised for SV-spanning reads.
  - Reference impl: `C++` · [philres/ngmlr](https://github.com/philres/ngmlr) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — same SW + chaining structure as minimap2
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-ngmlr`)
  - Consumes primitives: `block-aligner`, `noodles-bam`, future `rsomics-chain`
  - Notes: Pairs naturally with Sniffles (see [`sv-calling.md`](sv-calling.md)). Lower priority than minimap2 since modern Sniffles works fine on minimap2 output, but still used in SV-heavy pipelines. Code is small and approachable.

- [ ] **`Winnowmap`** — minimap2 derivative with weighted minimizers for repetitive references.
  - Reference impl: `C` · [marbl/Winnowmap](https://github.com/marbl/Winnowmap) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads (same as minimap2)
  - SIMD: inherits ksw2 SIMD from minimap2 fork
  - Quadrant: —
  - GPU-amenable: maybe — same kernel as minimap2
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-minimap` (small patch over minimap2 — natural fit as a preset)
  - Consumes primitives: same as minimap2 entry
  - Notes: Important for T2T-era human and plant genomes (centromeres, rDNA arrays). Mostly a small patch on top of minimap2; a Rust port should arrive once the minimap2 chainer is portable. Ship as a `--preset winnow` flag rather than a separate binary.

- [ ] **`LRA`** — long-read aligner with concurrent genomic + spliced modes.
  - Reference impl: `C++` · [ChaissonLab/LRA](https://github.com/ChaissonLab/LRA) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — same architectural opportunities as minimap2
  - Upstream license: `MIT`
  - Priority: `P2`
  - Layer: `B` (tool — `rsomics-lra` if user demand emerges)
  - Consumes primitives: `block-aligner`, `noodles-bam`
  - Notes: Niche; mostly Chaisson-lab structural variation work. Document interop only.

- [ ] **`lordFAST`** — PacBio CLR-focused aligner.
  - Reference impl: `C++` · [vpc-ccg/lordfast](https://github.com/vpc-ccg/lordfast) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — uncommon technology, low upside
  - Upstream license: `GPL-3.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: PacBio CLR is largely superseded by HiFi. GPL licence complicates derivative work. Skip unless a user demands it.
