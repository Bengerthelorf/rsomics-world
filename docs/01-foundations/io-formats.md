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
  - Reference impl: `C` (htslib) · [samtools/htslib](https://github.com/samtools/htslib) · `MIT/Expat`
  - Existing Rust: [`noodles`](https://github.com/zaeleus/noodles) (production, actively maintained)
  - Existing non-C alternatives: `htsjdk` (Java, used by GATK); `pysam` (Python wrapper of htslib)
  - Priority: `P0`
  - Notes: Adopt. Authoritative IO layer for the entire `rsomics-*` family.
    Contribute upstream rather than fork. Watch CRAM 3.1 codec compliance
    and async-tokio surface; both are improving but still flagged
    experimental.

- [x] **`needletail`** — fast FASTA/FASTQ parser with adapter detection.
  - Reference impl: `C` (kseq.h / readfq) · [lh3/readfq](https://github.com/lh3/readfq) · `MIT`
  - Existing Rust: [`needletail`](https://github.com/onecodex/needletail) (production)
  - Existing non-C alternatives: `seq_io` (Rust, lower-level)
  - Priority: `P0`
  - Notes: Adopt for FASTX hot paths. Internally uses `seq_io` algorithm.
    Faster than noodles-fastq on raw scan-only workloads because it skips
    UTF-8 validation. Use noodles when you need record-level metadata, use
    needletail when you need throughput.

- [x] **`rust-htslib`** — FFI bindings to htslib.
  - Reference impl: `C` · [samtools/htslib](https://github.com/samtools/htslib) · `MIT/Expat`
  - Existing Rust: [`rust-htslib`](https://github.com/rust-bio/rust-htslib) (`hts-sys` for raw bindings)
  - Existing non-C alternatives: `pysam`, `htsjdk`
  - Priority: `P1`
  - Notes: Transitional. Keep as a fallback when noodles lacks a feature
    (mpileup engine, some CRAM 3.1 corner cases). Each use site gets a
    tracking issue for migration to noodles.

- [~] **`SAM/BAM/CRAM`** — alignment record format family.
  - Reference impl: `C` · [samtools/htslib](https://github.com/samtools/htslib) · `MIT/Expat`
  - Existing Rust: `noodles-sam`, `noodles-bam`, `noodles-cram`
  - Existing non-C alternatives: `htsjdk` (Java)
  - Priority: `P0`
  - Notes: noodles covers SAM 1.6, BAM 1.6, CRAM 3.0 well. CRAM 3.1 codecs
    (rANS Nx16, fqzcomp) still maturing. Heavy users (Hi-C, single-cell)
    push noodles-bam edges first.

- [~] **`VCF/BCF`** — variant call format.
  - Reference impl: `C` · [samtools/bcftools](https://github.com/samtools/bcftools) · `MIT/Expat`
  - Existing Rust: `noodles-vcf`, `noodles-bcf`; also [`bcf_reader`](https://docs.rs/bcf_reader)
  - Existing non-C alternatives: `htsjdk` (Java); `cyvcf2` (Python/Cython)
  - Priority: `P0`
  - Notes: noodles is the path forward. Pop-gen scale (1000G, gnomAD) is a
    stress test: parsing 700M-variant BCFs needs careful streaming.

- [x] **`GFF/GTF`** — feature annotation formats.
  - Reference impl: `C++` · [The Sequence Ontology / Ensembl](http://gmod.org/wiki/GFF3) · `MIT-ish`
  - Existing Rust: [`noodles-gff`](https://crates.io/crates/noodles-gff),
    [`noodles-gtf`](https://crates.io/crates/noodles-gtf)
  - Existing non-C alternatives: `gffutils` (Python), `rtracklayer` (R)
  - Priority: `P0`
  - Notes: Adopt. Used by every quantifier and annotation tool.

- [x] **`BED`** — interval format.
  - Reference impl: `C++` · [arq5x/bedtools2](https://github.com/arq5x/bedtools2) · `MIT`
  - Existing Rust: `noodles-bed`; also `rust-bio` interval trees
  - Existing non-C alternatives: `pybedtools`
  - Priority: `P0`
  - Notes: Adopt noodles for parsing. The associated *operations* (intersect,
    merge, sort) belong to a future `rsomics-bedtools` crate — see module 09.

- [~] **`PAF`** — pairwise mapping format (minimap2 default output).
  - Reference impl: `C` · [lh3/minimap2](https://github.com/lh3/minimap2) · `MIT`
  - Existing Rust: [`paf`](https://crates.io/crates/paf) (parser only),
    `rustybam::paf`, `minimap2-paf-io`
  - Existing non-C alternatives: `paftools.js`
  - Priority: `P1`
  - Notes: Multiple small crates exist, none authoritative. A
    `noodles-paf` would resolve fragmentation. Track upstream PAF spec
    discussion under minimap2.

- [ ] **`MAF`** — multiple alignment format (UCSC) and Mutation Annotation
  Format (NCI/TCGA). Two *different* formats sharing a name; both unhandled.
  - Reference impl (UCSC): `C` · [UCSC kent tools](http://hgdownload.soe.ucsc.edu/admin/exe/) · `proprietary-but-free`
  - Reference impl (TCGA): `Python/Perl` · [mskcc/vcf2maf](https://github.com/mskcc/vcf2maf) · `Apache-2.0`
  - Existing Rust: none verified for either flavour
  - Existing non-C alternatives: `maftools` (R) for TCGA-MAF
  - Priority: `P2`
  - Notes: Niche but cited often in cancer pipelines. Start with TCGA-MAF
    because it interoperates with the VCF stack; UCSC-MAF can wait until
    comparative genomics work begins (module 08).

- [~] **`h5ad` / AnnData** — HDF5-backed single-cell matrix container.
  - Reference impl: `Python` · [scverse/anndata](https://github.com/scverse/anndata) · `BSD-3-Clause`
  - Existing Rust: [`anndata`](https://crates.io/crates/anndata) (subset),
    [`anndata-memory`](https://crates.io/crates/anndata-memory),
    [`single_rust`](https://crates.io/crates/single_rust),
    [`af-anndata`](https://crates.io/crates/af-anndata)
  - Existing non-C alternatives: `anndataR` (R); native HDF5 readers in Julia
  - Priority: `P0`
  - Notes: Fragmented landscape; no crate covers the full spec
    (`.layers`, `.raw`, `.uns` nested groups, `.obsm`/`.varm` arrays,
    backed/lazy access). Real target is a `rsomics-anndata` that the
    single-cell crates in module 04 all consume. Coordinate with SnapATAC2
    (uses anndata-rs internally) before forking.

- [ ] **`zarr`** — chunked array format (h5ad-zarr variant, spatial-omics
  next gen).
  - Reference impl: `Python` · [zarr-developers/zarr-python](https://github.com/zarr-developers/zarr-python) · `MIT`
  - Existing Rust: [`zarrs`](https://crates.io/crates/zarrs) (active, v3 spec)
  - Existing non-C alternatives: `zarr-java`, `tensorstore` (C++)
  - Priority: `P1`
  - Notes: Important for spatial transcriptomics + atlas-scale data. `zarrs`
    is promising; adopt if it stabilises before we hit module 04 spatial
    work.

- [ ] **`htsget`** — HTTP-streamed BAM/VCF.
  - Reference impl: `Java` · [ga4gh/htsget-refserver](https://github.com/ga4gh/htsget-refserver) · `Apache-2.0`
  - Existing Rust: `noodles-htsget` client; [`htsget-rs`](https://github.com/umccr/htsget-rs) server
  - Existing non-C alternatives: `htsget-refserver` (Java/Go)
  - Priority: `P2`
  - Notes: Cloud streaming. `htsget-rs` (UMCCR) is a strong Rust server
    implementation. Adopt; document interop with noodles client.
