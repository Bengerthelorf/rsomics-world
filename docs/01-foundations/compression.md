# Compression

> Codecs and CLI tools for the squashed bytes underneath every bioinformatics
> file.

## Scope

Block- and stream-level compression: DEFLATE/gzip, BGZF, zstd, lz4, xz, and
the user-facing CLIs (`bgzip`, `pigz`, `crabz`) that wrap them. Excludes
*record-level* file formats (those live in
[`io-formats.md`](io-formats.md)) and random-access indexes built on top of
BGZF (those live in [`indexing.md`](indexing.md)).

## Design notes

- BGZF is the workhorse: every BAM/VCF/CRAM/BCF/tabix file in production is
  a BGZF stream. Throughput here translates directly to pipeline wall time.
- Two implementation strategies coexist: pure-Rust DEFLATE via
  [`flate2`](https://github.com/rust-lang/flate2-rs) (with `miniz_oxide`
  or `zlib-ng` backends) and FFI to `libdeflate` via `libdeflater`. For
  block-sized inputs, libdeflate is ~2× faster — `gzp` and `crabz` use it
  by default.
- Multi-threaded compression is where Rust beats single-threaded `gzip` and
  matches `pigz`: see [`gzp`](https://github.com/sstadick/gzp) and
  [`crabz`](https://github.com/sstadick/crabz).
- zstd is a serious contender for *new* file formats (CRAM 3.1 uses it
  internally) but the existing bioinformatics ecosystem is overwhelmingly
  gzip/BGZF, so any zstd-only output needs a fallback path.
- xz/LZMA shows up only in archival contexts (e.g. SRA fasterq-dump
  outputs) and is not worth optimising.

## TODO

- [x] **`flate2`** — DEFLATE/gzip/zlib codec for Rust.
  - Reference impl: `C` · [madler/zlib](https://github.com/madler/zlib) · `Zlib`
  - Existing Rust: [`flate2`](https://github.com/rust-lang/flate2-rs)
    (pluggable: `miniz_oxide` pure-Rust, `zlib-ng` FFI, `cloudflare-zlib`)
  - Existing non-C alternatives: `sile/libflate` (pure-Rust DEFLATE)
  - Priority: `P0`
  - Notes: Adopt. The default `miniz_oxide` backend is portable and
    safe; switch to `zlib-ng` when throughput matters and ship a feature
    flag.

- [x] **`libdeflate`** — block-oriented DEFLATE optimised for known-size
  inputs.
  - Reference impl: `C` · [ebiggers/libdeflate](https://github.com/ebiggers/libdeflate) · `MIT`
  - Existing Rust: [`libdeflater`](https://crates.io/crates/libdeflater) (safe wrapper);
    [`libdeflate-sys`](https://crates.io/crates/libdeflate-sys) (raw)
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: FFI-only today. A pure-Rust libdeflate-equivalent (block-DEFLATE
    with vector intrinsics) is a real opportunity but a multi-month
    project. SIMD-critical inner loops.

- [x] **`BGZF`** — Blocked GNU Zip Format used by samtools/htslib.
  - Reference impl: `C` · [samtools/htslib/bgzf.c](https://github.com/samtools/htslib/blob/develop/bgzf.c) · `MIT/Expat`
  - Existing Rust: [`noodles-bgzf`](https://crates.io/crates/noodles-bgzf)
    (pure-Rust); [`bgzip`](https://crates.io/crates/bgzip);
    [`gzp`](https://crates.io/crates/gzp) (multithreaded)
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Adopt `noodles-bgzf` for IO correctness; pair with `gzp` /
    `libdeflater` for multithreaded *write* paths. noodles-bgzf reader is
    single-threaded for decompression — a parallel decoder is an open
    project (tracking issue
    [zaeleus/noodles#17](https://github.com/zaeleus/noodles/issues/17)).

- [x] **`bgzip` (CLI)** — `samtools` companion for creating BGZF files.
  - Reference impl: `C` · [samtools/htslib/bgzip.c](https://github.com/samtools/htslib) · `MIT/Expat`
  - Existing Rust: [`crabz`](https://github.com/sstadick/crabz) (pigz-style
    multithreaded gzip/BGZF CLI built on `gzp`)
  - Existing non-C alternatives: `bgzip` ships as part of htslib
  - Priority: `P1`
  - Notes: Adopt `crabz` as the `rsomics-bgzip` replacement. Already beats
    `bgzip --threads` on large inputs by using libdeflate per block.

- [x] **`pigz`** — parallel gzip CLI.
  - Reference impl: `C` · [madler/pigz](https://github.com/madler/pigz) · `Zlib`
  - Existing Rust: [`crabz`](https://github.com/sstadick/crabz);
    [`gzp`](https://crates.io/crates/gzp) (library)
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: `crabz` is the closest Rust equivalent. Ship as the default gzip
    binary in `rsomics-*` containers.

- [x] **`zstd`** — Facebook's zstandard codec.
  - Reference impl: `C` · [facebook/zstd](https://github.com/facebook/zstd) · `BSD-3-Clause / GPL-2.0`
  - Existing Rust: [`zstd`](https://github.com/gyscos/zstd-rs) (FFI),
    [`zstd-rs`](https://github.com/KillingSpark/zstd-rs) (pure-Rust decoder)
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Use FFI `zstd` for production (much faster, multi-threaded
    encoder). `KillingSpark/zstd-rs` is decoder-only and lags upstream; not
    yet a drop-in replacement.

- [x] **`lz4`** — fast streaming compressor.
  - Reference impl: `C` · [lz4/lz4](https://github.com/lz4/lz4) · `BSD-2-Clause`
  - Existing Rust: [`lz4_flex`](https://crates.io/crates/lz4_flex)
    (pure-Rust); `lz4-sys` (FFI)
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: `lz4_flex` is pure-Rust and within ~10% of the C library. Adopt.
    Used mostly for intermediate scratch files; rarely a user-facing format
    in genomics.

- [~] **`xz` / `liblzma`** — high-ratio LZMA codec.
  - Reference impl: `C` · [tukaani-project/xz](https://github.com/tukaani-project/xz) · `0BSD / LGPL`
  - Existing Rust: [`xz2`](https://crates.io/crates/xz2) (FFI);
    [`lzma-rs`](https://crates.io/crates/lzma-rs) (partial pure-Rust)
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: FFI is fine. Only relevant for SRA archive ingest and some
    legacy CRAM. Not a focus.

- [x] **`niffler`** — format-sniffing reader (auto-detect gzip/bgzf/zstd/xz).
  - Reference impl: — (Rust-native concept; similar to Python's `xopen`)
  - Existing Rust: [`niffler`](https://crates.io/crates/niffler);
    [`zopen`](https://crates.io/crates/zopen)
  - Existing non-C alternatives: `xopen` (Python)
  - Priority: `P1`
  - Notes: Adopt `niffler` as the default open-by-extension helper in
    CLI tools. Eliminates a class of "forgot to gunzip" user errors.
