# Transcript assembly and isoform reconstruction

> Reference-guided and de novo assembly of full-length transcripts /
> isoforms from RNA-seq reads.

## Scope

Reference-guided assemblers (StringTie, Cufflinks) start from a spliced
BAM and a genome; de novo assemblers (Trinity, Bridger, SOAPdenovo-Trans)
build transcripts directly from FASTQ without a reference. Long-read
isoform tools (IsoQuant) sit in between — they expect aligned long reads
but rely on graph traversal rather than read coverage models.

Quantification-only modes of these tools live in
[`quantification.md`](quantification.md).

## Design notes

- De novo transcript assembly is dominated by Trinity, which is
  Perl-orchestrated C++ — large, monolithic, hard to package. Two
  obvious points of leverage for Rust: (1) a clean orchestrator
  replacing `Trinity.pl` and (2) a faster Inchworm-equivalent k-mer
  extender.
- Reference-guided assembly is the more tractable rewrite. StringTie's
  flow network on the splice graph is well-defined and a good fit for
  `petgraph` + `noodles-bam`.
- Long-read isoform discovery is the future-facing target; IsoQuant
  is Python and CPU-bound, and a Rust rewrite could be 10x faster.
- We do not plan to revive Cufflinks; it is included only for legacy
  benchmark reproduction.

## TODO

- [ ] **`StringTie`** — flow-network reference-guided transcript
  assembler and quantifier.
  - Reference impl: `C++` · [gpertea/stringtie](https://github.com/gpertea/stringtie) · `MIT`
  - Existing Rust: none.
  - Existing non-C alternatives: none.
  - Priority: `P0`
  - Notes: Best target for a Rust assembler rewrite in this module.
    `petgraph` covers the splice graph data structure; flow algorithms
    (min-cost max-flow) need a small custom implementation.
    StringTie2 also handles long reads — implement both modes.

- [ ] **`Cufflinks`** — legacy reference-guided assembler.
  - Reference impl: `C++` · [cole-trapnell-lab/cufflinks](https://github.com/cole-trapnell-lab/cufflinks) · `Boost / OSS-friendly`
  - Existing Rust: none.
  - Existing non-C alternatives: StringTie is the official successor.
  - Priority: `P2`
  - Notes: Officially deprecated. List only so downstream Cuffdiff
    legacy pipelines have a documented replacement path. No porting
    intent.

- [ ] **`Trinity`** — de novo transcript assembly using the Inchworm /
  Chrysalis / Butterfly pipeline.
  - Reference impl: `C++ / Java / Perl` · [trinityrnaseq/trinityrnaseq](https://github.com/trinityrnaseq/trinityrnaseq) · `BSD-3-Clause`
  - Existing Rust: none.
  - Existing non-C alternatives: `rnaSPAdes` (C++, part of SPAdes).
  - Priority: `P1`
  - Notes: The whole pipeline is large; realistic Rust deliverables
    are (a) a `cargo`-installable orchestrator replacing the Perl
    driver, and (b) a Rust Inchworm using `rayon` + `kmer-rs` style
    rolling hashes. Don't try to rewrite Butterfly.

- [ ] **`Bridger`** — alternative de novo transcriptome assembler.
  - Reference impl: `C++` · [Bridger SourceForge mirror](https://sourceforge.net/projects/rnaseqassembly/) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: Trinity, rnaSPAdes.
  - Priority: `P2`
  - Notes: Cited as a Cufflinks-derived de novo assembler with full-length
    transcript reporting. Niche; not a priority rewrite.

- [ ] **`SOAPdenovo-Trans`** — BGI's de novo transcript assembler.
  - Reference impl: `C` · [aquaskyline/SOAPdenovo-Trans](https://github.com/aquaskyline/SOAPdenovo-Trans) · `GPL-3.0`
  - Existing Rust: none.
  - Existing non-C alternatives: Trinity.
  - Priority: `P2`
  - Notes: Useful for very large transcriptomes (plant / metazoan) due
    to lower memory footprint than Trinity. Niche.

- [ ] **`IsoQuant`** — long-read isoform discovery and quantification.
  - Reference impl: `Python` · [ablab/IsoQuant](https://github.com/ablab/IsoQuant) · `GPL-2.0`
  - Existing Rust: none.
  - Existing non-C alternatives: `bambu` (R), `FLAIR` (Python),
    `oarfish` (Rust; quantification only, no discovery).
  - Priority: `P1`
  - Notes: Pure-Python and slow on large nanopore datasets. A Rust
    rewrite that combines `minimap2-rs` + a splice-graph isoform
    discoverer would slot well with `oarfish` for quantification.
