# 01 — Foundations

Cross-cutting infrastructure every higher-level tool needs: file format I/O,
compression, indexing, building-block data structures, and the parallelism
patterns that hold it all together. Nothing in this module is biology-specific;
everything in modules 02–09 depends on this layer working well.

## Sub-docs

- [`io-formats.md`](io-formats.md) — FASTA/FASTQ, SAM/BAM/CRAM, VCF/BCF,
  GFF/GTF, BED, MAF, PAF, h5ad. Centred on `noodles` + `needletail`.
- [`compression.md`](compression.md) — gzip/bgzf/zstd/lz4/xz codecs and the
  `bgzip` / `pigz` CLI tools.
- [`indexing.md`](indexing.md) — fai/bai/csi/tbi/gzi random-access indexes
  and the `tabix` CLI.
- [`data-structures.md`](data-structures.md) — FM-index/BWT/suffix arrays,
  k-mer hashing (ntHash, MurmurHash3), MinHash, HyperLogLog, Bloom/Cuckoo
  filters.
- [`parallelism.md`](parallelism.md) — `rayon`, async I/O, GPU offload via
  `candle`/`burn`/`wgpu`; how big tools thread today and where Rust improves
  on that.

## Design posture

- Many of the canonical Rust crates here (`noodles`, `needletail`,
  `rust-htslib`) already exist and are production-grade. The work in this
  module is **audit, fill gaps, and document**, not "rewrite from zero".
- Where Rust *must* still ship a new crate (e.g. a pure-Rust bgzf writer
  that matches `libdeflate` throughput; an h5ad reader that does not depend
  on Python), the topic doc says so explicitly.
- All higher modules are forbidden from re-implementing format parsing.
  They consume the crates listed here.
