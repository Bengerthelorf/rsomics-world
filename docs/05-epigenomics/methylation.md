# DNA methylation

> Bisulfite-sequencing alignment / extraction and long-read modified-base
> calling.

## Scope

Two regimes:

1. **Short-read bisulfite sequencing (WGBS, RRBS)** — Bismark, BWA-meth,
   MethylDackel for extraction, methylpy / methylKit for analysis.
2. **Long-read modified-base calling (Nanopore, PacBio)** — modkit
   (modBAM-based) and the legacy nanopolish call-methylation.

Differential methylation testing at the regional level (DMR / DML) is
covered here in the analysis tools (methylKit, methylpy). Single-cell
methylation lives in module 04.

## Design notes

- The long-read side already has the canonical Rust tool: **`modkit`**
  from Oxford Nanopore. Adopt as-is.
- The short-read side is the most C / Perl-heavy corner of bioinformatics:
  Bismark is Perl wrapping Bowtie2, BWA-meth is a Python wrapper around
  BWA, MethylDackel is C++. There is a clear Rust opening for a single
  `rsomics-bisulfite` crate that subsumes alignment (via an in-silico
  C→T reference + a Rust short-read aligner from module 02) and
  extraction (per-read CpG / CHG / CHH counting).
- The `modBAM` (MM/ML tags) standard is now consolidated in the SAM
  specification and is supported by `noodles-bam` — modkit reads/writes
  it natively. A future bisulfite extractor should emit modBAM too,
  for unified downstream tooling.
- `methylpy` and `methylKit` are downstream analysis layers (DMR finding,
  smoothing). methylKit is R + C++; methylpy is Python. Neither is
  algorithmically large; both can be wrapped via `extendr` / pyo3 and
  the few hot loops (smoothing, beta-binomial likelihood) ported into
  rsomics for free speedups.

## TODO

- [ ] **`Bismark`** — Perl-orchestrated bisulfite aligner over Bowtie2.
  - Reference impl: `Perl` · [FelixKrueger/Bismark](https://github.com/FelixKrueger/Bismark) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: BWA-meth (Python+BWA)
  - Parallelism: upstream Perl wrapping Bowtie2 pthreads
  - SIMD: inherited from bowtie2
  - Quadrant: —
  - GPU-amenable: maybe — same SW kernel rationale as short-read aligners
  - Upstream license: `GPL-3.0`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-bisulfite` (alignment mode flag)
  - Consumes primitives: `rsomics-fm-index`, `block-aligner`, `noodles-bam`, `noodles-fasta`
  - Notes: Pure-Rust short-read bisulfite aligner is a credible target once `02-genomics` ships a Rust short-read aligner. Strategy: in-silico convert both reference and reads (C→T forward, G→A reverse), align with `rsomics-align`, restore base identities and emit a Bismark-compatible BAM.

- [ ] **`BWA-meth`** — Python wrapper around BWA for bisulfite alignment.
  - Reference impl: `Python + C` · [brentp/bwa-meth](https://github.com/brentp/bwa-meth) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream BWA pthreads
  - SIMD: inherited from BWA
  - Quadrant: —
  - GPU-amenable: maybe — same as Bismark
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-bisulfite` (default alignment mode)
  - Consumes primitives: `rsomics-bwa`, `block-aligner`, `noodles-bam`, `noodles-fasta`
  - Notes: Lighter than Bismark and faster in benchmarks. Rust port rides on the same `rsomics-align` foundation. Brent Pedersen's other tools (`echtvar`, `vcfexpress`) are already Rust — natural upstream contact.

- [ ] **`methylKit`** — R/Bioconductor methylation analysis (DMR / DML).
  - Reference impl: `R / C++` · [Bioconductor methylKit](https://bioconductor.org/packages/release/bioc/html/methylKit.html) · `Artistic-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: R BiocParallel + C++ inner loops
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — DMR testing is per-region parallel
  - Upstream license: `Artistic-2.0`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-bisulfite` (DMR/DML mode + extendr bridge)
  - Consumes primitives: `extendr`-bridge, `noodles-bed`, future `rsomics-stats` (logistic-regression, smoothing)
  - Notes: Bridge via `extendr`. The few inner loops (logistic-regression DMR testing, smoothing) port well to `rsomics-bisulfite` if benchmarks show a bottleneck.

- [ ] **`methylpy`** — Python methylation extractor + DMR finder.
  - Reference impl: `Python / C++` · [yupenghe/methylpy](https://github.com/yupenghe/methylpy) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: no — Python orchestration over per-region work
  - Upstream license: `MIT`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-bisulfite`
  - Consumes primitives: `noodles-bam`, future `rsomics-stats`
  - Notes: Niche; mostly used in Ecker-lab pipelines. Wrap via PyO3 if needed; rewrite is low priority.

- [ ] **`MethylDackel`** — universal methylation extractor for bisulfite BAMs.
  - Reference impl: `C` · [dpryan79/MethylDackel](https://github.com/dpryan79/MethylDackel) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: no — per-position binomial counting, memory-latency-bound
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `subcommand-of-rsomics-bisulfite` (extraction mode — primary user-facing subcommand)
  - Consumes primitives: `noodles-bam`, `noodles-bed`, `statrs` (binomial)
  - Notes: Small (< 5 kLoC), focused, MIT-licensed, and called by every nf-core methylseq pipeline. Excellent Rust port target — `noodles-bam` covers IO, the per-position binomial logic is tiny. Output format (CpG bedGraph + methylKit) is well-specified.

- [x] **`modkit`** — Oxford Nanopore modified-base toolkit (modBAM-based).
  - Reference impl: `Rust` · [nanoporetech/modkit](https://github.com/nanoporetech/modkit) · `MPL-2.0`
  - Existing Rust: [`modkit`](https://github.com/nanoporetech/modkit) (binary tool, install from source — not on crates.io under that name; `cf-modkit-macros` on crates.io is unrelated)
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: —
  - Parallelism: rayon
  - SIMD: inherits htslib's hand SIMD via `rust-htslib`
  - Quadrant: ②
  - GPU-amenable: no — bedMethyl summarisation is memory-bandwidth-bound
  - Upstream license: `MPL-2.0`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Adopt as-is. Canonical tool for nanopore modBAM → bedMethyl and DMR/DML analysis on long reads. The modkit *algorithm* is rust-native (ONT's own modBAM processing), but the BAM IO depends on `rust-htslib` (FFI wrapper of C htslib), putting the effective perf class at Quadrant ② — same kind-vs-quadrant distinction as Sawfish and oarfish. The effective quadrant becomes ① if/when modkit migrates to `noodles-bam`. crates.io name not published — install from source (squat catalog updated).

- [ ] **`nanopolish`** (call-methylation) — legacy HMM modified-base caller from signal-level data.
  - Reference impl: `C++` · [jts/nanopolish](https://github.com/jts/nanopolish) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: modkit (modBAM-based, supersedes signal-level calling for current basecallers)
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: no — HMM signal-level, legacy
  - Upstream license: `MIT`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Mostly superseded by Dorado + modkit. Keep only for legacy R9 nanopore data analysis. Not a porting target.
