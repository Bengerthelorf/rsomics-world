# 06 — Metagenomics

> Classification, profiling, MAG-resolved assembly, and amplicon analysis of
> shotgun and 16S/ITS sequencing data from microbial communities.

## Sub-areas

| File | Scope |
|------|-------|
| [`classification.md`](classification.md) | Read-level taxonomic classification: Kraken2/KrakenUniq, Centrifuge, Kaiju, MetaPhlAn4, mOTUs, sourmash, Ganon, MetaMaps, CCMetagen. |
| [`profiling.md`](profiling.md) | Quantitative community + functional profiling: HUMAnN3, Bracken, MicrobeCensus, MetaPhlAn4 (profiling face), eggNOG-mapper, KEGG/GhostKOALA. |
| [`assembly-mag.md`](assembly-mag.md) | Metagenome assembly, binning, bin refinement, dereplication, QC, and GTDB taxonomy: MEGAHIT, metaSPAdes, IDBA-UD, MetaBAT2, MaxBin2, CONCOCT, SemiBin2, VAMB, DAS_Tool, CheckM/CheckM2, GTDB-Tk, dRep. |
| [`amplicon.md`](amplicon.md) | 16S/ITS amplicon pipelines: QIIME2, DADA2, mothur, VSEARCH/USEARCH, swarm, Deblur, PICRUSt2. |

## Design notes

- Metagenomics is dominated by k-mer indexing (Kraken2-style) and sketch
  comparisons (sourmash, MinHash). Rust already has a flagship implementation
  here — `sourmash`'s core is Rust — so this is the most "Rust-ready" of all
  the specialty domains.
- Classifier index sizes are the dominant constraint, not CPU. Pure-Rust
  rewrites need to attack memory layout and mmap discipline, not just
  algorithmic tricks. `noodles` + memory-mapped index crates (`memmap2`)
  give us the right primitives.
- Binning is shifting from feature-based clustering (TNF + coverage) to deep
  learning (SemiBin2, VAMB, TaxVamb). Rust's path forward is through `candle`
  or `burn` for inference; training will stay on PyTorch for the near term.
- HUMAnN3 is essentially a glue pipeline around DIAMOND + MetaPhlAn4 — the
  Rust win here is replacing the Python orchestration layer, not the
  alignment kernels.
- Amplicon analysis is fragmented across Python (QIIME2), R (DADA2, phyloseq),
  and standalone C (VSEARCH, swarm). A pure-Rust DADA2-style denoiser would
  be a high-value, scoped target with a clean math core.
- License watch: Kraken2, MEGAHIT, Bracken, GTDB-Tk, MetaBAT2, raxml-ng, and
  many adjacent tools are **GPL**. Any derivative or close port must respect
  that — `CONVENTIONS.md` requires flagging GPL inheritance in the relevant
  module doc; we re-state it in each topic file.
