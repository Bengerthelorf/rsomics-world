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
  - Existing Rust: none.
  - Existing non-C alternatives: BWA-meth (Python+BWA), bwa-meth-rs has
    not appeared.
  - Priority: `P1`
  - Notes: Pure-Rust short-read bisulfite aligner is a credible target
    once `02-genomics` ships a Rust short-read aligner. Strategy:
    in-silico convert both reference and reads (C→T forward, G→A
    reverse), align with `rsomics-align`, restore base identities and
    emit a Bismark-compatible BAM.

- [ ] **`BWA-meth`** — Python wrapper around BWA for bisulfite alignment.
  - Reference impl: `Python + C` · [brentp/bwa-meth](https://github.com/brentp/bwa-meth) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Lighter than Bismark and faster in benchmarks. Rust port
    rides on the same `rsomics-align` foundation. Brent Pedersen's
    other tools (`echtvar`, `vcfexpress`) are already Rust — natural
    upstream contact.

- [ ] **`methylKit`** — R/Bioconductor methylation analysis (DMR / DML).
  - Reference impl: `R / C++` · [Bioconductor methylKit](https://bioconductor.org/packages/release/bioc/html/methylKit.html) · `Artistic-2.0`
  - Existing Rust: none.
  - Priority: `P1`
  - Notes: Bridge via `extendr`. The few inner loops (logistic-regression
    DMR testing, smoothing) port well to `rsomics-methyl` if benchmarks
    show a bottleneck.

- [ ] **`methylpy`** — Python methylation extractor + DMR finder.
  - Reference impl: `Python / C++` · [yupenghe/methylpy](https://github.com/yupenghe/methylpy) · `MIT`
  - Existing Rust: none.
  - Priority: `P2`
  - Notes: Niche; mostly used in Ecker-lab pipelines. Wrap via PyO3 if
    needed; rewrite is low priority.

- [ ] **`MethylDackel`** — universal methylation extractor for
  bisulfite BAMs.
  - Reference impl: `C` · [dpryan79/MethylDackel](https://github.com/dpryan79/MethylDackel) · `MIT`
  - Existing Rust: none.
  - Priority: `P0`
  - Notes: Small (< 5 kLoC), focused, MIT-licensed, and called by every
    nf-core methylseq pipeline. Excellent Rust port target —
    `noodles-bam` covers IO, the per-position binomial logic is tiny.
    Output format (CpG bedGraph + methylKit) is well-specified.

- [x] **`modkit`** — Oxford Nanopore modified-base toolkit (modBAM-based).
  - Reference impl: `Rust` · [nanoporetech/modkit](https://github.com/nanoporetech/modkit) · `MPL-2.0`
  - Existing Rust: modkit (this row).
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Adopt as-is. Canonical tool for nanopore modBAM →
    bedMethyl and DMR/DML analysis on long reads. Already uses
    `noodles-bam`; rsomics packages it and contributes when needed.

- [ ] **`nanopolish`** (call-methylation) — legacy HMM modified-base
  caller from signal-level data.
  - Reference impl: `C++` · [jts/nanopolish](https://github.com/jts/nanopolish) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: modkit (modBAM-based, supersedes
    signal-level calling for current basecallers).
  - Priority: `P2`
  - Notes: Mostly superseded by Dorado + modkit. Keep only for legacy
    R9 nanopore data analysis. Not a porting target.
