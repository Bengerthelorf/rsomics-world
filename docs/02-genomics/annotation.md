# Variant annotation

> Functional, frequency, and curated annotation of called variants.

## Scope

Tools that *add information* to a VCF: consequence prediction (VEP, SnpEff,
ANNOVAR, Funcotator), population-frequency / database joins (echtvar,
slivar), filtering / expression languages over annotated VCFs. Excludes
the upstream variant calling (see [`variant-calling.md`](variant-calling.md))
and downstream prioritisation / clinical reporting (a future module).

## Design notes

- The big three (VEP, SnpEff, ANNOVAR) all attempt the same primary task
  — predicting consequence of variants against transcript models — and
  produce *concordant only ~65%* of the time. Choice of annotator
  materially changes downstream filtering. Any Rust rewrite must be
  benchmarked against the originals on a fixed gnomAD-class input.
- The bottleneck in modern annotation pipelines is **database joining**,
  not consequence prediction. echtvar's million-variant-per-second
  performance shows that Rust + compact integer encoding is the right
  approach. slivar uses the same insight (genotype-aware filtering at
  scale).
- VEP's transcript model + Ensembl integration is the hardest part to
  replicate. A native Rust GFF/GTF-driven consequence engine is
  feasible but requires careful handling of edge cases (NMD, frameshift
  rescue, MANE select, refseq vs Ensembl coordinates).
- For new work, the lightest-weight target is a `rsomics-snpeff`-style
  consequence predictor — algorithm is well-understood, output format
  (VCF ANN field) is standardised, and `noodles-{vcf,gff,fasta}` cover
  most of the IO.

## TODO

- [ ] **`VEP` (Variant Effect Predictor)** — Ensembl's consequence
  annotator.
  - Reference impl: `Perl` · [Ensembl/ensembl-vep](https://github.com/Ensembl/ensembl-vep) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Slowest of the big three but most feature-rich (plugins,
    custom annotations, cache architecture, Ensembl-MANE integration).
    A Rust reimplementation should aim at the **consequence prediction
    core** first, leaving plugins and the Perl ecosystem alone. Cache
    format compatibility is a major engineering concern.

- [ ] **`SnpEff`** — Java consequence annotator.
  - Reference impl: `Java` · [pcingola/SnpEff](https://github.com/pcingola/SnpEff) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Simpler internal model than VEP, more amenable to clean-room
    re-derivation. MIT licence. Output ANN field is the *de-facto*
    annotation standard for many community pipelines.

- [ ] **`ANNOVAR`** — table-based annotator.
  - Reference impl: `Perl` · [annovar.openbioinformatics.org](https://annovar.openbioinformatics.org/) · proprietary academic licence
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Licence is restrictive (academic-only, no redistribution). We
    cannot port it. Document interop only: read its table-based outputs
    via `noodles`/`polars`.

- [ ] **`Funcotator`** — GATK's annotator.
  - Reference impl: `Java` · [broadinstitute/gatk](https://github.com/broadinstitute/gatk) (Funcotator subpackage) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Mostly used inside GATK-flavoured pipelines. Outputs VCF or
    MAF; the MAF writer is uniquely useful (TCGA workflows). Lower
    priority than VEP/SnpEff but ships MAF support.

- [x] **`echtvar`** — fast compressed annotation joiner.
  - Reference impl: `Rust` · [brentp/echtvar](https://github.com/brentp/echtvar) · `MIT`
  - Existing Rust: [`echtvar`](https://github.com/brentp/echtvar)
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Adopt. Million-variant-per-second annotation lookups using
    compact integer encodings; benchmark winner for population-frequency
    joins. Already production in many gnomAD pipelines.

- [ ] **`slivar`** — VCF filtering + trio analysis with embedded JS.
  - Reference impl: `Nim` · [brentp/slivar](https://github.com/brentp/slivar) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Nim-based, not C. Brentp pairs slivar + echtvar in his own
    pipelines. A Rust port (`rsomics-slivar`?) could share the echtvar
    backing store and offer an expression DSL via `rhai` / `rune`. Worth
    consolidating with `vcfexpress` (next entry).

- [x] **`vcfexpress`** — expression language for VCF filtering (Rust).
  - Reference impl: `Rust` · [brentp/vcfexpress](https://github.com/brentp/vcfexpress) · `MIT`
  - Existing Rust: [`vcfexpress`](https://github.com/brentp/vcfexpress)
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Adopt. Lua-embedded VCF filter language; pairs with echtvar
    and `noodles-vcf`. Underdocumented but the right architectural model
    for fast filterable annotation pipelines.

- [ ] **`bcftools csq`** — bcftools' consequence subcommand.
  - Reference impl: `C` · [samtools/bcftools](https://github.com/samtools/bcftools) · `MIT/Expat`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Lightweight, haplotype-aware (handles MNVs and adjacent
    variants together). The cleanest reference implementation in the
    consequence-prediction space. Worth porting as part of the
    `rsomics-bcftools` umbrella.
