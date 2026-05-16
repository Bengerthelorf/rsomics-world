# rsomics-seqio

Fast FASTQ input reader for the rsomics-* tool family. Layer A foundation crate.

## Problem

`needletail` is hardwired to flate2 (miniz_oxide) and decompresses gzip on the
same thread that processes records, so decode and parse never overlap and the
decompressor is the slowest available. fastp instead uses a dedicated reader
thread with large block buffers and Intel ISA-L igzip. `rsomics-seqio` matches
that architecture and goes one step further: the reader thread does *only*
decompression, and record parsing runs in parallel across the rayon pool.

## Design

`open_fastq(path)` detects the input format by magic bytes (never by file
extension) and dispatches to one of three paths:

1. **Plain** — `BufReader<File>` + line parser. No extra threads.
2. **Plain gzip** — a dedicated producer thread decompresses in 8 MiB blocks
   and sends whole-record-aligned raw byte slabs (`Vec<u8>`) over a bounded
   `crossbeam-channel`. The consumer parses each slab in parallel on the rayon
   thread pool, so per-record allocation is off the reader thread. The producer
   never parses; the consumer never decompresses.
3. **BGZF** (feature `bgzf`) — `noodles-bgzf` parallel block decoder wrapped in
   the same line parser.

Slab boundary invariant: every slab the producer emits contains an integer
number of complete 4-line FASTQ records. Boundaries are found by counting
newlines in groups of four, which is immune to `@`/`+` appearing as quality
characters. The tail after the last complete record is carried into the next
block. Truncated or corrupt input surfaces as `Err`, never a short clean read;
a panic in the decode thread is delivered to the consumer as a loud `Err`.

## Decompressor backends (gz path)

| Feature | Backend | Quadrant | Notes |
|---|---|---|---|
| `igzip-backend` (default) | `rsomics-igzip` (Intel ISA-L igzip via `isal-sys`) | ② FFI over C+asm | fastp-class 4 MiB-in / 8 MiB-out block shape; requires nasm + a C toolchain |
| `pure` | flate2 + zlib-rs | ① pure Rust | Degraded fallback; no nasm/C toolchain; NOT for production publishing |
| `bgzf` | noodles-bgzf 0.47 | ① pure Rust | BGZF only; noodles parallel block decode |

All `unsafe` lives in the isolated `rsomics-igzip` crate (Quadrant ②); this
crate stays 100% safe and keeps `[lints] workspace = true`. The producer /
consumer code is identical across backends — only the `build_decoder` call
differs.

## BGZF detection

Detection uses the gzip EXTRA-field `BC` subfield per SAM spec §4.1, not the
`0x1f 0x8b` magic alone or the file extension. A plain-gz file lacking the BC
subfield is classified as Gz, never as Bgzf — misclassification would corrupt
reads. The header probe fills across short reads so a partial first `read`
cannot misclassify the input.

## OwnedRecord

```rust
pub struct OwnedRecord {
    pub id:   Vec<u8>,   // header after '@', line-terminator stripped, otherwise verbatim
    pub seq:  Vec<u8>,   // sequence line, line-terminator stripped, otherwise verbatim
    pub qual: Vec<u8>,   // quality line, line-terminator stripped, otherwise verbatim
}
```

`id` semantics match needletail's `SequenceRecord::id()` (full header after `@`,
trailing `\n` and a preceding `\r` stripped, whitespace not split), so a tool
migrating from needletail to `rsomics_seqio::open_fastq` sees zero output
divergence. The parser enforces `seq.len() == qual.len()`.

## Origin

This crate is an independent Rust implementation; no upstream source was
copied. The reader-thread + large-block buffer architecture is documented in
the fastp paper (Chen et al. 2018, doi:10.1093/bioinformatics/bty560) and the
fastp source (MIT, `src/fastqreader.{h,cpp}`, `src/peprocessor.cpp`). The
parallel-parse split is original.

The default gz backend is Intel ISA-L igzip, wrapped behind the isolated
`rsomics-igzip` crate (the only place in the workspace where audited `unsafe`
lives). BGZF format: SAM/BAM specification §4.1
(https://samtools.github.io/hts-specs/SAMv1.pdf).

Upstream credits:
- fastp (MIT) — reader-thread + large-block architecture reference
- Intel ISA-L (BSD-2-Clause), via `rsomics-igzip` / `isal-sys` — default gz backend
- noodles-bgzf (MIT) — BGZF decode

License: MIT OR Apache-2.0.
