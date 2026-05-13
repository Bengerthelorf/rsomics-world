# Metagenomic classification

> Read-level taxonomic assignment from shotgun metagenomic sequencing,
> covering exact-match k-mer classifiers, FM-index / BWT classifiers,
> protein-space classifiers, marker-gene tools, and sketch-based methods.

## Scope

Includes: per-read or per-contig taxonomic label assignment, both nucleotide
(Kraken2-style) and protein (Kaiju, MMseqs2-style); marker-gene profilers
that double as classifiers (MetaPhlAn4, mOTUs); long-read aware classifiers
(MetaMaps, sourmash on PacBio/ONT). Excludes: abundance re-estimation
(see [profiling](profiling.md) for Bracken, HUMAnN3) and assembly-time
binning (see [assembly-mag](assembly-mag.md) for MetaBAT2 etc.).

## Design notes

- Kraken2's database is a perfect-hash map from minimizer-space k-mers to LCA
  taxa. A pure-Rust rewrite is mostly a careful `fxhash`/`ahash` + mmap +
  zero-copy parser exercise; the algorithm itself is well documented.
- Centrifuge and Kaiju both ride on FM-indices. The `fm-index` crate exists
  but is not at production scale; substantial work needed before we replace
  them. Worth checking `bio-types` and `rust-bio`'s suffix-array tooling.
- MetaPhlAn4 and mOTUs are marker-gene tools — most of their value is the
  curated marker DB and the strain-resolution math, not raw CPU. A Rust
  rewrite is feasible but the leverage is in the database curation pipeline.
- `sourmash` is already pure-Rust at the core with Python bindings. We adopt
  unconditionally and consider contributing rather than forking.
- Ganon2 uses Interleaved Bloom Filters — a data structure Rust handles
  cleanly with bit-packed arrays; this is a tractable rewrite target.
- License watch: Kraken2 family (Kraken2/KrakenUniq/Bracken) is **GPL-3**;
  Kaiju is **GPL-3**; mOTUs is **GPL-3**; MetaPhlAn is MIT. Clean-room
  pure-Rust ports avoid the derivative-work question; direct ports inherit
  the GPL.

## TODO

- [ ] **`kraken2`** — exact-match minimizer k-mer classifier; the de-facto baseline.
  - Reference impl: `C++` · [DerrickWood/kraken2](https://github.com/DerrickWood/kraken2) · `MIT`
  - Existing Rust: none mature; some experimental `kraken2-rs` repos exist but none verified production-ready
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `ganon2` (C++ with IBF, more modern data structure)
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — minimizer hashing is SIMT-trivial, hash lookup latency-bound
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-kraken`)
  - Consumes primitives: `rsomics-kmer` (nthash-rs minimizer rolling hash), `noodles-fastq`, `memmap2`, `fxhash`/`ahash`, `rayon`
  - Notes: One of the most-cited bioinformatics tools of the past decade. Memory-mapped k-mer→taxon hash is the core; Rust's `memmap2` + `fxhash` + SIMD-friendly minimizer rolling hash (`nthash-rs` crate) give us everything we need. Build `rsomics-kraken` that can both read upstream `.k2d` indexes and emit its own format.

- [ ] **`KrakenUniq`** — Kraken extended with HyperLogLog unique-k-mer counting for specificity.
  - Reference impl: `C++` · [fbreitwieser/krakenuniq](https://github.com/fbreitwieser/krakenuniq) · `MIT / GPL-3` (dual; Kraken-1-derived parts GPL-3)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — same as kraken2 plus HLL register merge
  - Upstream license: `MIT / GPL-3` (dual)
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-kraken` (a `--unique-counts` mode)
  - Consumes primitives: `rsomics-kraken`, `amadeus-streaming` (HLL)
  - Notes: Once `rsomics-kraken` exists, add a HyperLogLog sketch per taxon. Rust HLL crates are mature (`amadeus-streaming` has explicit SIMD). Low marginal cost once the base classifier is done.

- [ ] **`Centrifuge`** — BWT/FM-index metagenomic classifier from Salzberg group.
  - Reference impl: `C++` · [DaehwanKimLab/centrifuge](https://github.com/DaehwanKimLab/centrifuge) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `Centrifuger` (C++ successor by the same group, more memory-efficient)
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — FM-index probing same as short-read aligners
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: `B` (tool — `rsomics-centrifuge`)
  - Consumes primitives: `rsomics-fm-index`, `noodles-fastq`, future `rsomics-stats`
  - Notes: Largely superseded by `Centrifuger` and by k-mer methods. Lower priority than Kraken2. If we port, target `Centrifuger`'s data structures rather than legacy Centrifuge.

- [ ] **`MetaPhlAn4`** — clade-specific marker-gene profiler.
  - Reference impl: `Python` (wraps Bowtie2) · [biobakery/MetaPhlAn](https://github.com/biobakery/MetaPhlAn) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing + Bowtie2 pthreads
  - SIMD: inherits Bowtie2
  - Quadrant: —
  - GPU-amenable: maybe — marker alignment is SW-like
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-metaphlan`)
  - Consumes primitives: `minimap2` or future `rsomics-bowtie`, `noodles-bam`, `polars`, future `rsomics-stats` (strain-level)
  - Notes: Hot loop is Bowtie2 alignment against ~1M marker DB. Replace with `minimap2-rs` short-read mode or a custom SeedExtend on marker-DB index. Python orchestration is the big rewrite win. Strain-level (StrainPhlAn) is a separate, harder follow-up.

- [ ] **`mOTUs`** (motus 3/4) — marker-gene OTU profiler.
  - Reference impl: `Python` (wraps BWA-MEM) · [motu-tool/mOTUs](https://github.com/motu-tool/mOTUs) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing + BWA pthreads
  - SIMD: inherits BWA
  - Quadrant: —
  - GPU-amenable: maybe — alignment SIMT-friendly
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-metaphlan` (same marker-gene umbrella; --tool motus flag)
  - Consumes primitives: `rsomics-bwa`, `noodles-bam`, `polars`
  - Notes: Same shape as MetaPhlAn4 — Python orchestration around an aligner + curated DB. Worth doing only if the aligner crate in `rsomics-align` becomes the standard, then mOTUs and MetaPhlAn4 become thin wrappers.

- [x] **`sourmash`** — FracMinHash sketching for sequence comparison and gather.
  - Reference impl: `Rust` core + `Python` CLI · [sourmash-bio/sourmash](https://github.com/sourmash-bio/sourmash) · `BSD-3-Clause`
  - Existing Rust: [`sourmash`](https://crates.io/crates/sourmash) `0.22.0` (the official core)
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: `branchwater` (Rust multithreaded plugin from the same group)
  - Parallelism: rayon
  - SIMD: auto-vectorize (underlying hashes vectorise)
  - Quadrant: ①
  - GPU-amenable: maybe — sketch comparison parallelises
  - Upstream license: `BSD-3-Clause`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Production Rust. Use it directly from any downstream crate that needs sketch-based search/gather/cluster. Consider upstream contributions if we hit limits.

- [ ] **`Kaiju`** — protein-space (BWT) metagenomic classifier.
  - Reference impl: `C++` · [bioinformatics-centre/kaiju](https://github.com/bioinformatics-centre/kaiju) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `MMseqs2` (C++, broader feature set; protein search is a strict superset of Kaiju)
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — protein-space BWT probing
  - Upstream license: `GPL-3`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-mmseqs` (a `--mode kaiju` flag within the protein-search umbrella)
  - Consumes primitives: `rsomics-fm-index`, future `rsomics-mmseqs`, `noodles-fastq`
  - Notes: Protein-space classification matters for highly diverged or novel organisms. The right Rust target is probably a focused MMseqs2-prefilter port rather than reimplementing Kaiju's full BWT.

- [ ] **`CCMetagen`** — KMA-alignment based eukaryote/prokaryote classifier.
  - Reference impl: `Python` (wraps KMA, which is C) · [vrmarcelino/CCMetagen](https://github.com/vrmarcelino/CCMetagen) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python + KMA pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: no — small user base, no upside
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Niche but valuable for fungal/eukaryotic detection. Bottleneck is KMA, not Python glue. Skip unless we already have a KMA-equivalent crate.

- [ ] **`Ganon` / `Ganon2`** — Interleaved Bloom Filter classifier with built-in DB management.
  - Reference impl: `C++` + `Python` glue · [pirovc/ganon](https://github.com/pirovc/ganon) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — IBF lookup is SIMT-trivial
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-ganon`)
  - Consumes primitives: `rsomics-kmer`, `fastbloom`, `roaring` (bit-vector), `polars` for DB metadata
  - Notes: IBF is one of the cleaner data structures to implement in Rust; `roaring`/bit-vector crates exist. The DB-management story (incremental updates of NCBI taxonomy) is half the value. Good scoped target.

- [ ] **`MetaMaps`** — long-read approximate mapping + EM composition estimation.
  - Reference impl: `C++` · [DiltheyLab/MetaMaps](https://github.com/DiltheyLab/MetaMaps) · `Public Domain`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `sourmash gather` for long reads; Centrifuger long-read mode
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — same as minimap2-style approximate mapping
  - Upstream license: `Public Domain`
  - Priority: `P2`
  - Layer: `subcommand-of-sourmash` (or contribute long-read composition EM to sourmash/branchwater upstream)
  - Consumes primitives: `sourmash`, `branchwater`, `minimap2-rs`
  - Notes: Built on MashMap-style approximate mapping. With `minimap2-rs` and `sourmash` available, an entirely new MetaMaps port is hard to justify; better to extend `sourmash`/`branchwater` for long-read composition EM.
