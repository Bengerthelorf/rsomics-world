# Design principles

Rules of thumb every crate under the `rsomics-*` umbrella should follow.

## 1. Streaming first, batch second

Bioinformatics inputs are big. Default APIs should iterate, not load. A
`fasta::Reader` returns an `Iterator<Item = Record>`, not a `Vec<Record>`. If
a user wants random access, they ask for it (via an index).

## 2. Zero-copy where it does not hurt readability

Records borrow from the underlying buffer when possible (`&[u8]` over
`Vec<u8>`). Allocation-per-record loops are the most common reason naive
Rust ports lose to C — be deliberate.

## 3. SIMD is opt-in but discoverable

Inner loops (Smith-Waterman, hash lookup, base counting) get
SIMD-accelerated paths gated behind `cfg(target_feature)` or runtime
dispatch. The scalar path is always present and tested. We do **not** ship
crates that require AVX-512 to build.

## 4. `rayon` for data parallelism, `tokio` rarely

Most bioinformatics workloads are CPU-bound, embarrassingly parallel chunks
of records. `rayon::par_iter()` is the default tool. `tokio` and `async`
appear only where I/O genuinely dominates (cloud storage clients, web
servers).

## 5. No global state, no panics in libraries

Library crates return `Result<_, E>` for everything fallible. `panic!` is
for invariants the caller cannot violate. Global mutable state is banned
(no `static mut`, no `lazy_static` registries).

## 6. Match the data, not the GUI

Output formats are the standard ones: BAM, VCF, BED, GFF, MTX, h5ad/zarr.
We do not invent file formats. If a tool's original output is messy, we
match it bit-for-bit by default and offer a `--clean` flag for the
sanitized variant.

## 7. Wrap before rewriting

If a mature C library exists and its license allows, we may ship a
`*-sys` FFI wrapper crate alongside the pure-Rust port during the
transition. The wrapper is *temporary scaffolding* — every wrapper has a
tracking issue for replacement.

## 8. Benchmark or it didn't happen

Every crate ships a `benches/` directory with `criterion` benchmarks
against (a) the C original via FFI or subprocess and (b) a reference
dataset (1000 Genomes chr22, PBMC 10k, etc.). PRs that change perf-critical
code show the delta.

## 9. CLI and library are the same crate, or sibling crates

`rsomics-vcf` is the library. `rsomics-bcftools` is the CLI built on top
(if needed). The library never `println!`s or `std::process::exit`s.

## 10. Document for the biologist, not the compiler

`cargo doc` output is the user manual. Examples in docstrings use real
datasets where possible. Error messages name the file and offset, not just
"invalid input."
