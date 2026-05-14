# rsomics-bam

CLI over the [`rust-htslib`][rust-htslib] FFI bindings to [`htslib`][htslib]
for BAM / SAM / CRAM operations. **FFI wrapper, not a clean-room rewrite** —
the actual binary I/O is htslib C code; this crate provides a uniform Rust
CLI surface compatible with `samtools`'s expected behavior.

Quadrant ② in the workspace convention (`docs/00-overview/`): we adopt the
upstream C library through Rust bindings rather than re-implementing it.
htslib is the perf and correctness reference for BAM I/O; we don't compete
with it on those axes.

## Install

```sh
cargo install rsomics-bam
```

System prerequisites for the build: `libz` (zlib) and a C toolchain are
required by the `hts-sys` crate's C build. `bzip2`, `lzma`, and `curl`
are NOT required at v0.1.0 (the corresponding rust-htslib features are
disabled).

## Usage

```sh
# Count records in a BAM (analogous to samtools view -c)
rsomics-bam view -c input.bam

# Full --help
rsomics-bam --help
rsomics-bam view --help
```

Subcommands shipped in 0.1.0:

| subcommand | description | analog |
|---|---|---|
| `view -c` | count records in a BAM | `samtools view -c` |

`view` without `--count`, plus `sort`, `index`, `markdup`, region
filtering, SAM output, and CRAM support land in later 0.1.x / 0.2.x
releases.

## Compatibility

The `view -c` count is verified against `samtools view -c` on every CI
push (Linux runners). See `tests/view_test.rs` for the assertion.

The crate uses rust-htslib's high-level API; semantics (record-count
definition, header skip rules, EOF detection) match htslib exactly
because htslib is what executes the read.

## Origin

This crate is an **FFI wrapper** over the public C library
[`htslib`][htslib] via the Rust bindings crate
[`rust-htslib`][rust-htslib]. It is not a clean-room reimplementation —
htslib's C code does the actual BAM/SAM/CRAM parsing, compression, and
indexing.

- Upstream htslib: <https://github.com/samtools/htslib> (MIT)
- Rust bindings: <https://github.com/rust-bio/rust-htslib> (MIT)
- Reference behaviour: `samtools` <https://github.com/samtools/samtools> (MIT)

License: MIT OR Apache-2.0 (this crate). Upstream htslib + rust-htslib are
MIT; the FFI relationship inherits that compatibly.

The motivation for adopting rather than re-implementing: htslib is over
twenty years of carefully-optimised C with extensive bug history and
correctness audits. A pure-Rust rewrite would be a multi-year project
(BGZF, BAM index, CRAM, MD/NM/cigar handling, every flag edge case) and
would race against htslib's perf floor for no clear ecosystem benefit
at this scale. Future Quadrant ① pure-Rust ports may emerge for specific
hot paths (e.g., BGZF decompression via `std::simd`); they will live in
foundation crates and be opt-in.

[rust-htslib]: https://github.com/rust-bio/rust-htslib
[htslib]: https://github.com/samtools/htslib
