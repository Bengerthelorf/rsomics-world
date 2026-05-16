# rsomics-igzip

Minimal Quadrant-② FFI wrapper over Intel ISA-L igzip for fast streaming gzip decompression.

Exposes a single safe public type: `GzReader`, which implements `std::io::Read`.
All `unsafe` is contained within this crate.  Every consumer (`rsomics-seqio`,
`rsomics-fastq-*`) stays 100% safe Rust and keeps `[lints] workspace = true`.

## Why this crate exists

| Crate | Problem |
|---|---|
| `isal-rs` | Safe wrapper but hard-codes `BUF_SIZE = 16 KiB`, preventing the large-block read pattern that gives ISA-L its throughput advantage |

The correct shape — 4 MiB compressed input / 8 MiB decompressed output blocks,
mirroring fastp `src/fastqreader.cpp readToBufIgzip` with multi-member gzip
support — requires calling `isal_sys` directly with caller-controlled buffer
sizes.

## Platform support

ISA-L's hand-written aarch64 assembly does not assemble under Apple's
integrated assembler, so `isal-sys` is a Linux-only dependency. On non-Linux
targets this crate still compiles but `GzReader::new` returns an `Unsupported`
error; consumers (e.g. `rsomics-seqio`) select a pure-Rust decoder per target.
The performance contract is enforced on Linux.

## Dependency quadrant

| Crate | Quadrant | Why |
|---|---|---|
| `rsomics-igzip` | ② FFI | `build.rs` compiles ISA-L C + NASM via `cc`; isal-sys is the raw `-sys` crate |

## Origin

Independent minimal Rust FFI wrapper over Intel ISA-L igzip via `isal-sys`.

- **ISA-L**: Intel Intelligent Storage Acceleration Library,
  <https://github.com/intel/isa-l>, BSD-2-Clause.
- **isal-sys**: <https://github.com/milesgranger/isal-rs>, MIT.
- **Algorithmic reference shape**: fastp `src/fastqreader.cpp` `readToBufIgzip`,
  S. Chen et al., *iMeta* 2023, <https://github.com/OpenGene/fastp>, MIT.
  Read-permitted under MIT; cited as the algorithmic reference, not copied verbatim.

No GPL source was used.

## License

MIT OR Apache-2.0
