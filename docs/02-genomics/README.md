# 02 — Genomics

DNA-centric workflows: read alignment, assembly, variant calling (small +
structural), variant annotation, and preprocessing. The bread-and-butter of
short-read and long-read genomics pipelines.

## Sub-docs

- [`alignment-short-read.md`](alignment-short-read.md) — BWA, Bowtie2, SNAP,
  Strobealign, NovoAlign.
- [`alignment-long-read.md`](alignment-long-read.md) — minimap2, NGMLR,
  Winnowmap, LRA, lordFAST.
- [`assembly.md`](assembly.md) — SPAdes, ABySS, MaSuRCA (short-read);
  Flye, hifiasm, Canu, wtdbg2, NextDenovo, Shasta, Raven, Verkko
  (long-read); polishers (Pilon, Racon, Medaka, NextPolish).
- [`variant-calling.md`](variant-calling.md) — GATK HaplotypeCaller,
  DeepVariant, FreeBayes, bcftools call, Strelka2, Octopus, Clair3,
  PEPPER-Margin-DeepVariant, DeepTrio.
- [`sv-calling.md`](sv-calling.md) — Manta, Delly, Lumpy, SvABA, GRIDSS
  (short-read); Sniffles, cuteSV, SVIM, pbsv, Sawfish, Wham (long-read).
- [`annotation.md`](annotation.md) — VEP, SnpEff, ANNOVAR, slivar, echtvar,
  Funcotator.
- [`preprocessing.md`](preprocessing.md) — fastp, Trimmomatic, cutadapt,
  BBDuk, Trim Galore, AdapterRemoval, FastQC, MultiQC, seqkit, seqtk.

## Design posture

- Highest density of mature C/C++ tools in the entire stack. Rust ports
  are mostly *aspirational* here — only a few production-grade examples
  exist (varlociraptor for variant calling; echtvar for annotation;
  GGCAT for cDBG assembly). Most other tools require new pure-Rust crates
  built on the foundations layer.
- SIMD/SWAR matters enormously: BWA-MEM2's main improvement over BWA-MEM
  is AVX-512 chaining; minimap2's inner loop is hand-vectorised. Any
  rewrite needs a credible SIMD plan.
- We treat **`bwa-mem2`** (not original `bwa`) as the modern-C++ baseline
  for short-read alignment performance comparisons. The 2-3× speedup over
  bwa-mem is the bar to clear in pure Rust.
- For long-read alignment we likely *adopt* `minimap2-rs` (FFI) and only
  attempt a pure-Rust port later — the algorithm is moving target and the
  cost/benefit is poor.
- For variant calling, the deep-learning callers (DeepVariant, Clair3,
  PEPPER) need `candle`/`burn` integration, not a from-scratch model
  training run. We ship inference, not training.
