# I/O formats

> Parsers, writers, and record types for the file formats every other
> module reads or writes.

## Scope

Covers the textual and binary file formats that move sequence, alignment,
variant, feature, and single-cell data between tools: FASTA/FASTQ,
SAM/BAM/CRAM, VCF/BCF, GFF/GTF, BED, MAF, PAF, h5ad. The boundary with the
neighbouring topic [`compression.md`](compression.md) is: this doc covers
*record-level* readers/writers; the codecs themselves (DEFLATE, BGZF, zstd)
live one level down. Random-access indexes (fai/bai/csi) live in
[`indexing.md`](indexing.md).

## Design notes

- Rust is in a strong position here. [`noodles`](https://github.com/zaeleus/noodles)
  is a pure-Rust, spec-tracking implementation of nearly every sequencing
  format and is the de-facto reference library. Most "rewrites" in this
  topic are not new crates but contributions back to noodles.
- The single largest hole is **h5ad / AnnData**. Several crates (`anndata`,
  `anndata-memory`, `single_rust`/`SingleRust`) exist but each implements
  only a subset of the spec; a unified pure-Rust AnnData layer is still
  needed for single-cell workflows.
- Streaming-first APIs matter: a typical aligned BAM is 100+ GB. Every
  reader in this layer must expose `Iterator<Item = Record>` rather than
  forcing `Vec<Record>`.
- Zero-copy record views over the input buffer (`&[u8]` fields, not
  `String` everywhere) are the main differentiator vs. naive ports of
  htslib's API.
- For tools that already wrap htslib (`rust-htslib`), treat them as
  transitional. Goal is to make `noodles` strictly superset their feature
  coverage so we can drop the FFI dependency.

## TODO

- [x] **`noodles`** — pure-Rust SAM/BAM/CRAM/VCF/BCF/GFF/GTF/BED/FASTA/FASTQ/BGZF/CSI/tabix.
  - Reference impl: `C` · [samtools/htslib](https://github.com/samtools/htslib) · `MIT`
  - Existing Rust: [`noodles`](https://github.com/zaeleus/noodles) `0.110.0` (workspace meta-crate; subcrates carry their own versions)
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: `htsjdk` (Java, GATK backbone); `pysam` (CPython binding of htslib)
  - Parallelism: streaming-iterator (`Iterator<Item = Record>`); rayon-amenable on the consumer side
  - SIMD: auto-vectorize at record-decode level; explicit SIMD comes via codec deps (`zlib-rs`, `libdeflater`)
  - Quadrant: ①
  - GPU-amenable: no — parsing is I/O-bound, not compute
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Authoritative IO layer for the entire `rsomics-*` family. Contribute upstream rather than fork. Watch CRAM 3.1 codec compliance and async-tokio surface; both are improving but still flagged experimental. Edition 2024, MSRV 1.89, workspace under one repo.

- [x] **`needletail`** — fast FASTA/FASTQ parser with adapter detection.
  - Reference impl: `C` · [lh3/readfq](https://github.com/lh3/readfq) (kseq.h) · `MIT`
  - Existing Rust: [`needletail`](https://github.com/onecodex/needletail) `0.7.3`
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: `seq_io` (Rust, lower-level)
  - Parallelism: streaming-iterator; designed for outer-loop rayon parallelism per chunk
  - SIMD: auto-vectorize; relies on bytewise scan (`memchr`) rather than handwritten SIMD
  - Quadrant: ①
  - GPU-amenable: no — sequential parsing, I/O-bound
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Adopt for FASTX hot paths. Internally uses `seq_io` algorithm. Faster than `noodles-fastq` on raw scan-only workloads because it skips UTF-8 validation. Use noodles when you need record-level metadata, needletail when you need throughput.

- [x] **`rust-htslib`** — FFI bindings to htslib.
  - Reference impl: `C` · [samtools/htslib](https://github.com/samtools/htslib) · `MIT`
  - Existing Rust: [`rust-htslib`](https://github.com/rust-bio/rust-htslib) `1.0.0` (paired with [`hts-sys`](https://crates.io/crates/hts-sys) `2.2.0` for raw bindings)
  - Existing Rust kind: `FFI-wrapper`
  - Existing non-C alternatives: `pysam`, `htsjdk`
  - Parallelism: inherits htslib's pthread model on `bam_mt_*` paths; otherwise single-threaded
  - SIMD: inherits htslib's compile-time SIMD (CRC32C / popcnt intrinsics)
  - Quadrant: ②
  - GPU-amenable: no — htslib is CPU-only
  - Upstream license: `MIT` (htslib also `MIT`)
  - Priority: `P1`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Transitional. Keep as a fallback when noodles lacks a feature (mpileup engine, some CRAM 3.1 corner cases). Each use site gets a tracking issue for migration to noodles. The 1.0 release in 2026-04 makes it a stable transitional target.

- [~] **`SAM/BAM/CRAM`** — alignment record format family.
  - Reference impl: `C` · [samtools/htslib](https://github.com/samtools/htslib) · `MIT`
  - Existing Rust: [`noodles-sam`](https://crates.io/crates/noodles-sam) `0.85.0`, [`noodles-bam`](https://crates.io/crates/noodles-bam) `0.89.0`, [`noodles-cram`](https://crates.io/crates/noodles-cram) `0.93.0`
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: `htsjdk` (Java)
  - Parallelism: streaming-iterator + per-record async tokio variants
  - SIMD: auto-vectorize; BGZF codec deps carry explicit SIMD
  - Quadrant: ①
  - GPU-amenable: no — record-level parsing, not compute
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: noodles covers SAM 1.6 and BAM 1.6 in full; CRAM 3.1 codecs (rANS Nx16, fqzcomp) are landing but still flagged "experimental" in some entry points. Heavy users (Hi-C, single-cell) push noodles-bam edges first.

- [~] **`VCF/BCF`** — variant call format.
  - Reference impl: `C` · [samtools/bcftools](https://github.com/samtools/bcftools) · `MIT`
  - Existing Rust: [`noodles-vcf`](https://crates.io/crates/noodles-vcf) `0.88.0`, [`noodles-bcf`](https://crates.io/crates/noodles-bcf) `0.86.0`; supplementary [`bcf_reader`](https://github.com/bguo068/bcf-reader) `0.3.2`
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: `htsjdk` (Java); `cyvcf2` (Python/Cython)
  - Parallelism: streaming-iterator for noodles; `bcf_reader` adds explicit rayon over records
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: no — record-level parsing, not compute
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: noodles is the default; `bcf_reader` is a lightweight option for pop-gen scale (1000G, gnomAD) where parsing 700M-variant BCFs needs careful streaming and the rayon-over-records pattern wins.

- [x] **`GFF/GTF`** — feature annotation formats.
  - Reference impl: `C++` · [The Sequence Ontology / Ensembl](http://gmod.org/wiki/GFF3) · spec is public domain; reference parsers in `gffread` (`MIT`)
  - Existing Rust: [`noodles-gff`](https://crates.io/crates/noodles-gff) `0.57.0`, [`noodles-gtf`](https://crates.io/crates/noodles-gtf) `0.52.0`
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: `gffutils` (Python), `rtracklayer` (R)
  - Parallelism: streaming-iterator
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: no — text parsing, I/O-bound
  - Upstream license: `MIT` (for gffread; spec itself is open)
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Used by every quantifier and annotation tool downstream.

- [x] **`BED`** — interval format.
  - Reference impl: `C++` · [arq5x/bedtools2](https://github.com/arq5x/bedtools2) · `MIT`
  - Existing Rust: [`noodles-bed`](https://crates.io/crates/noodles-bed) `0.34.0`; supplementary `rust-bio` interval trees
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: `pybedtools`
  - Parallelism: streaming-iterator; the *operations* (intersect/merge/sort) are rayon-amenable when implemented
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: no — text parsing; operations themselves are interval-tree work, CPU-friendly
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `adopt` (parsing); the operations belong to `rsomics-bedtools` (Layer B, module 09)
  - Consumes primitives: —
  - Notes: Adopt noodles for parsing. Operations crate is downstream of `rsomics-intervals` (foundation, [`data-structures.md`](data-structures.md)).

- [~] **`PAF`** — pairwise mapping format (minimap2 default output).
  - Reference impl: `C` · [lh3/minimap2](https://github.com/lh3/minimap2) · `MIT`
  - Existing Rust: [`paf`](https://github.com/ARU-life-sciences/paf) `0.2.1` (parser only, borderline stale 2024-10); `rustybam::paf` (part of a binary toolkit)
  - Existing Rust kind: `partial-port`
  - Existing non-C alternatives: `paftools.js` (JS distributed with minimap2)
  - Parallelism: single-threaded in the existing parsers
  - SIMD: none
  - Quadrant: ④
  - GPU-amenable: no — small text parser
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `A` (foundation — a future `rsomics-paf` or contribution to `noodles-paf`)
  - Consumes primitives: —
  - Notes: `paf` crate last update 2024-10-29 is borderline stale — kept because PAF spec is stable and parsers don't need active dev. `rustybam` is a binary tool, not a clean library. Fragmented landscape; a `noodles-paf` would resolve it. Track upstream PAF spec discussion under minimap2.

- [ ] **`MAF`** — multiple alignment format (UCSC) and Mutation Annotation Format (NCI/TCGA). Two *different* formats sharing a name; both unhandled.
  - Reference impl (UCSC): `C` · [UCSC kent tools](http://hgdownload.soe.ucsc.edu/admin/exe/) · UCSC academic source license (free for non-commercial)
  - Reference impl (TCGA): `Python/Perl` · [mskcc/vcf2maf](https://github.com/mskcc/vcf2maf) · `Apache-2.0`
  - Existing Rust: none verified for either flavour
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `maftools` (R) for TCGA-MAF
  - Parallelism: single-threaded in upstreams
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — text parsing
  - Upstream license: see above (two different flavours)
  - Priority: `P2`
  - Layer: `A` (foundation — future `rsomics-maf`)
  - Consumes primitives: —
  - Notes: Niche but cited often in cancer pipelines. Start with TCGA-MAF because it interoperates with the VCF stack; UCSC-MAF can wait until comparative genomics work begins (module 08).

- [~] **`h5ad` / AnnData** — HDF5-backed single-cell matrix container.
  - Reference impl: `Python` · [scverse/anndata](https://github.com/scverse/anndata) · `BSD-3-Clause`
  - Existing Rust: [`anndata`](https://github.com/kaizhang/anndata-rs) `0.6.2` (workspace; subset of spec); [`anndata-memory`](https://github.com/SingleRust/Anndata-Memory) `1.0.7`; [`single_rust`](https://github.com/SingleRust/SingleRust) `0.5.8`; [`af-anndata`](https://github.com/COMBINE-lab/af-anndata) `0.4.1`
  - Existing Rust kind: `partial-port` (each covers a different subset)
  - Existing non-C alternatives: `anndataR` (R); native HDF5 readers in Julia
  - Parallelism: `anndata-rs` and `anndata-memory` lean on `ndarray + rayon`; `single_rust` also uses `nalgebra + rayon` and optional tokio
  - SIMD: auto-vectorize via ndarray; HDF5 codec layer (blosc) carries SIMD
  - Quadrant: ①+② (in-memory ops are ①, HDF5 IO is ② via `hdf5-metno-sys` + `blosc-src` + `libz-sys`)
  - GPU-amenable: maybe — only at the on-device array layer (Vectorize / Candle); HDF5 IO is FFI-bound and stays on CPU
  - Upstream license: `BSD-3-Clause`
  - Priority: `P0`
  - Layer: `A` (foundation — future `rsomics-anndata`)
  - Consumes primitives: future `rsomics-stats` (compute side); HDF5 IO is FFI-bound today
  - Notes: Fragmented landscape; no crate covers the full spec (`.layers`, `.raw`, `.uns` nested groups, `.obsm`/`.varm` arrays, backed/lazy access). Real target is `rsomics-anndata` that the single-cell crates in module 04 all consume. Coordinate with SnapATAC2 (uses anndata-rs internally) before forking. The HDF5 FFI dependency is the architectural challenge — a pure-Rust HDF5 reader for the AnnData subset would let us cross into Quadrant ① end-to-end.

- [ ] **`zarr`** — chunked array format (h5ad-zarr variant, spatial-omics next gen).
  - Reference impl: `Python` · [zarr-developers/zarr-python](https://github.com/zarr-developers/zarr-python) · `MIT`
  - Existing Rust: [`zarrs`](https://github.com/zarrs/zarrs) `0.23.11` (active, v3 spec)
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: `zarr-java`, `tensorstore` (C++)
  - Parallelism: explicit rayon (`rayon`, `rayon_iter_concurrent_limit`)
  - SIMD: explicit (optional `simd-adler32` checksum codec)
  - Quadrant: ①
  - GPU-amenable: maybe — at the on-device array layer; not at the chunked IO layer itself
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: `zarrs` is the strongest pure-Rust Zarr v3 implementation; build.rs is a thin metadata generator (no C compile). Important for spatial transcriptomics and atlas-scale data. Adopt when module 04 spatial work begins.

- [ ] **`htsget`** — HTTP-streamed BAM/VCF.
  - Reference impl: `Java` · [ga4gh/htsget-refserver](https://github.com/ga4gh/htsget-refserver) · `Apache-2.0`
  - Existing Rust: [`noodles-htsget`](https://crates.io/crates/noodles-htsget) `0.11.0` (client); [`htsget-rs`](https://github.com/umccr/htsget-rs) (server, edition 2024 workspace with actix/axum/lambda backends)
  - Existing Rust kind: `rust-native` (independent Rust impls of the GA4GH htsget spec; Java refserver is one of several spec implementations, not a C/C++ upstream to port)
  - Existing non-C alternatives: `htsget-refserver` (Java/Go)
  - Parallelism: tokio async for both client and server
  - SIMD: inherits codec SIMD via BGZF deps
  - Quadrant: ①
  - GPU-amenable: no — HTTP streaming
  - Upstream license: `Apache-2.0` (spec); `MIT` (htsget-rs)
  - Priority: `P2`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: `htsget-rs` (UMCCR) is a strong Rust server implementation. Adopt; document interop with the noodles client.
