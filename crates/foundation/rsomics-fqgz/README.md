# rsomics-fqgz

Layer-A FASTQ output writer shared by the rsomics-* `fastq-*` tools.

`ChunkedWriter` writes FASTQ records to a file. When the path ends in `.gz`
it emits one self-contained gzip member per ~256 KB chunk, compressing the
chunks in parallel with libdeflate; otherwise it writes plain text.

## Why

gzip permits concatenated members, and `gunzip` / fastp / seqkit decode them
transparently — so per-chunk members can be compressed independently across
the rayon pool. libdeflate is the fast compressor here: flate2/zlib-rs is
~3× slower single-threaded than fastp's libdeflate, so a pure-Rust deflate
would be the bottleneck on the write side.

This is the single home for that machinery: it was previously copied into
each `fastq-*` tool; a module needed by 2+ tool crates belongs in Layer A.

## Origin

Independent Rust implementation. libdeflate via the `libdeflater` crate
(Quadrant ②, FFI over the libdeflate C library, MIT). The per-chunk
concatenated-gzip-member layout matches fastp's writer behaviour (fastp,
MIT, Chen et al. 2018, doi:10.1093/bioinformatics/bty560).

License: MIT OR Apache-2.0.
