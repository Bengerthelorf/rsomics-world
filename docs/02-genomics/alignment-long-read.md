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

- [~] **`minimap2`** — versatile pairwise aligner, the *de-facto* long-read
  aligner.
  - Reference impl: `C` · [lh3/minimap2](https://github.com/lh3/minimap2) · `MIT`
  - Existing Rust: [`minimap2-rs`](https://github.com/jguhlin/minimap2-rs) (FFI wrapper);
    [`minimap2-temp`](https://crates.io/crates/minimap2-temp);
    no pure-Rust port verified
  - Existing non-C alternatives: `mm2-fast` (C++/SIMD fork by bwa-mem2 group)
  - Priority: `P0`
  - Notes: FFI-only Rust today. Pure-Rust port is the long-term goal but
    is a multi-person-year effort: chainer, ksw2-like SW, presets, splice
    model. Adopt the FFI wrapper for now; track `mm2-fast` upstream
    improvements.

- [ ] **`NGMLR`** — long-read mapper specialised for SV-spanning reads.
  - Reference impl: `C++` · [philres/ngmlr](https://github.com/philres/ngmlr) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Pairs naturally with Sniffles (see
    [`sv-calling.md`](sv-calling.md)). Lower priority than minimap2 since
    modern Sniffles works fine on minimap2 output, but still used in
    SV-heavy pipelines. Code is small and approachable.

- [ ] **`Winnowmap`** — minimap2 derivative with weighted minimizers for
  repetitive references.
  - Reference impl: `C` · [marbl/Winnowmap](https://github.com/marbl/Winnowmap) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Important for T2T-era human and plant genomes (centromeres,
    rDNA arrays). Mostly a small patch on top of minimap2; a Rust
    port should arrive once the minimap2 chainer is portable.

- [ ] **`LRA`** — long-read aligner with concurrent genomic + spliced modes.
  - Reference impl: `C++` · [ChaissonLab/LRA](https://github.com/ChaissonLab/LRA) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Niche; mostly Chaisson-lab structural variation work.
    Document interop only.

- [ ] **`lordFAST`** — PacBio CLR-focused aligner.
  - Reference impl: `C++` · [vpc-ccg/lordfast](https://github.com/vpc-ccg/lordfast) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: PacBio CLR is largely superseded by HiFi. GPL licence
    complicates derivative work. Skip unless a user demands it.
