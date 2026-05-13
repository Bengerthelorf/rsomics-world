# Transcript and gene quantification

> Turning aligned (or pseudo-aligned) reads into transcript / gene abundance
> estimates.

## Scope

Two families:

1. **Lightweight / selective-alignment quantifiers** (Salmon, kallisto,
   piscem) — index a transcriptome and produce TPM/counts directly from
   FASTQ.
2. **Counters on top of BAM** (featureCounts, HTSeq-count, RSEM) — start
   from a STAR/HISAT2 BAM and a GTF.

StringTie does both quantification and assembly; the assembly half lives
in [`assembly-isoform.md`](assembly-isoform.md).

## Design notes

- Selective alignment + EM is the area where Rust already has the most
  visible presence: COMBINE-lab's `oarfish` is pure Rust, `simpleaf` is
  pure Rust, and `piscem` is a Rust binary wrapping a C++ index. Adopting
  these and contributing back is more productive than greenfield porting.
- For BAM-based counters the bottleneck is interval-tree querying of GFF
  features. `rust-bio` has interval trees, `noodles-bam` has efficient
  iteration, so a `featureCounts` clone is a small project with high
  pipeline impact.
- HTSeq-count is Python and slow but still cited because of strict-mode
  semantics. A drop-in Rust replacement with bit-identical counts in
  intersection-strict / intersection-nonempty / union modes would be
  immediately useful.
- RSEM does EM over a transcriptome-aligned BAM. Its math is shared with
  Salmon / oarfish, but its IO layer is dated and a Rust rewrite is
  attractive.

## TODO

- [ ] **`Salmon`** — selective-alignment-based transcript quantification.
  - Reference impl: `C++` · [COMBINE-lab/salmon](https://github.com/COMBINE-lab/salmon) · `GPL-3.0`
  - Existing Rust: none directly; the successor `piscem` is a Rust binary around a C++ core, see below. `alevin-fry` (Rust) consumes salmon / piscem RAD output.
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `piscem` (Rust+C++ hybrid)
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — EM iteration is dense linear algebra; selective-alignment scoring is SW-like
  - Upstream license: `GPL-3.0`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-salmon`)
  - Consumes primitives: `rsomics-kmer`, `rsomics-fm-index`, `block-aligner`, `noodles-fastq`, `noodles-fasta`, `noodles-gff`, `noodles-bam`, future `rsomics-stats` (EM)
  - Notes: Bulk-quant feature parity with `salmon quant` is the single most impactful Rust deliverable in transcriptomics. The maths (selective-alignment scoring + range-factor EM) is well-documented; the engineering pain is the compact suffix-array index. Coordinate with COMBINE-lab on whether `piscem` becomes the canonical path.

- [~] **`kallisto`** — pseudoalignment-based transcript quantification using the T-DBG.
  - Reference impl: `C++` · [pachterlab/kallisto](https://github.com/pachterlab/kallisto) · `BSD-2-Clause`
  - Existing Rust: partial — [`debruijn_mapping`](https://github.com/10XGenomics/rust-pseudoaligner) `0.6.0` (10x's repo; binary tool, install from source — not on crates.io despite the package name in its Cargo.toml) implements the T-DBG pseudoalignment primitive in Rust, used inside Cell Ranger but not packaged as a kallisto drop-in
  - Existing Rust kind: `partial-port`
  - Existing non-C alternatives: `BUStools` (companion, C++)
  - Parallelism: rayon
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: maybe — T-DBG traversal is irregular, EM is GPU-friendly
  - Upstream license: `BSD-2-Clause`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-kallisto`)
  - Consumes primitives: `debruijn` (`rust-debruijn` ecosystem), `rsomics-kmer`, `noodles-fastq`, future `rsomics-stats` (EM), future `rsomics-bus` for BUS-format IO
  - Notes: `debruijn_mapping` plus an EM-quantification layer would give most of `kallisto quant`. Output formats (BUS, h5) need `noodles` / custom HDF5 layers (`hdf5-metno`). Both names `rust-pseudoaligner` (the repo) and `debruijn_mapping` (the Cargo package name) are **unpublished on crates.io** — install from `10XGenomics/rust-pseudoaligner` source.

- [ ] **`RSEM`** — EM-based transcript expression estimation from transcriptome alignments.
  - Reference impl: `C++ / Perl` · [deweylab/RSEM](https://github.com/deweylab/RSEM) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `salmon --alignment-mode` covers most use cases
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — EM iteration is dense linear algebra
  - Upstream license: `GPL-3.0`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-rsem`)
  - Consumes primitives: `noodles-bam`, `noodles-fasta`, future `rsomics-stats` (EM), `ndarray`
  - Notes: Still mandatory for some isoform-DE pipelines that expect `.isoforms.results`. A Rust port would be a relatively small project once `noodles` + `nalgebra`/`ndarray` are in place. Match RSEM's log-likelihood model exactly so DESeq2 / tximport downstream is unaffected.

- [ ] **`featureCounts`** (Subread) — interval-based read summarization over a GFF.
  - Reference impl: `C` · [Subread/featureCounts](https://subread.sourceforge.net/featureCounts.html) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `htseq-count` (Python)
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: no — interval-tree querying is memory-latency-bound
  - Upstream license: `GPL-3.0`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-featurecounts`)
  - Consumes primitives: `noodles-bam`, `noodles-gff`, `rsomics-intervals` (foundation)
  - Notes: Most commonly-used gene-level counter in bulk RNA-seq. Small, focused tool — implement as a `noodles-bam` + interval-tree binary. Match output column-for-column with `Rsubread::featureCounts` so existing DESeq2 / edgeR scripts work unchanged.

- [ ] **`HTSeq-count`** — Python read-feature counter, reference for intersection-strict semantics.
  - Reference impl: `Python / Cython` · [htseq/htseq](https://github.com/htseq/htseq) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `featureCounts` (faster, different semantics)
  - Parallelism: Python serial / multiprocessing
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — same constraint as featureCounts
  - Upstream license: `GPL-3.0`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-featurecounts` (a `--mode htseq-strict` flag inside the same binary)
  - Consumes primitives: `noodles-bam`, `noodles-gff`, `rsomics-intervals`
  - Notes: Slow; pipelines tolerate it for strict-mode semantics. A Rust drop-in matching `--mode intersection-strict` + `--mode union` bit-for-bit is a useful Phase-2 deliverable. Folds into `rsomics-featurecounts` as a counting-semantics flag.

- [ ] **`StringTie`** (quantification mode) — transcript abundance over spliced BAMs.
  - Reference impl: `C++` · [gpertea/stringtie](https://github.com/gpertea/stringtie) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — splice-graph EM is dense linear algebra
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-stringtie` (same binary covers assembly + quant; see [`assembly-isoform.md`](assembly-isoform.md))
  - Consumes primitives: `noodles-bam`, `noodles-gff`, `rsomics-intervals`, future `rsomics-stats`
  - Notes: See `assembly-isoform.md` for the assembly side; the `-e` (estimate only, given GTF) mode is just quantification. Same rewrite covers both.

- [x] **`oarfish`** — long-read transcript quantification with coverage-aware EM.
  - Reference impl: `Rust` · [COMBINE-lab/oarfish](https://github.com/COMBINE-lab/oarfish) · `BSD-3-Clause`
  - Existing Rust: [`oarfish`](https://crates.io/crates/oarfish) `0.9.2`
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: `IsoQuant` (Python), `bambu` (R)
  - Parallelism: rayon
  - SIMD: auto-vectorize on EM iteration
  - Quadrant: ② (uses `minimap2-rs` for alignment in its hot path; the EM layer is pure-Rust ①, but the alignment dep is FFI)
  - GPU-amenable: maybe — EM iteration is dense linear algebra; alignment is the minimap2 question
  - Upstream license: `BSD-3-Clause`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Adopt as-is. COMBINE-lab uses `minimap2-rs` underneath. We package it in the rsomics ecosystem and contribute upstream when needed. The Quadrant ② reflects the FFI alignment dep; when the future pure-Rust minimap2 ports land, oarfish's effective quadrant becomes ①.

- [x] **`piscem`** — next-gen pure-Rust selective-alignment + index for Salmon / alevin-fry.
  - Reference impl: `Rust` (originally `Rust + C++` upstream; the C++ core has been replaced by pure-Rust `piscem-rs` + `cf1-rs`) · [COMBINE-lab/piscem](https://github.com/COMBINE-lab/piscem) · `BSD-3-Clause`
  - Existing Rust: [`piscem`](https://crates.io/crates/piscem) `0.20.0` (Rust CLI front-end); [`piscem-rs`](https://crates.io/crates/piscem-rs) `0.5.0` (mapping engine, pure-Rust reimplementation of the original C++ piscem); [`cf1-rs`](https://crates.io/crates/cf1-rs) `0.4.0` (compacted reference dBG index, pure-Rust)
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: salmon, kallisto
  - Parallelism: rayon throughout
  - SIMD: auto-vectorize; no FFI codec dep
  - Quadrant: ①
  - GPU-amenable: maybe — selective-alignment scoring is SW-like, EM is dense
  - Upstream license: `BSD-3-Clause`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Adopt. The piscem 0.20.0 era is fully pure-Rust — `piscem-rs` is the pure-Rust reimplementation of the original C++ mapping engine, `cf1-rs` is the pure-Rust compacted reference-dBG index. The `build.rs` in the front-end binary only reads `Cargo.lock` to embed `cf1-rs`'s version into the output; it does **not** compile C++. This is the closest existing analogue to a future `rsomics-salmon` and the strongest pure-Rust quantification stack today.
