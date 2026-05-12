# Motivation

## Why bioinformatics needs new tooling

Most workhorse bioinformatics software was written between 2008 and 2014, in C
or C++, by graduate students or small academic groups. They are extraordinary
achievements (BWA, samtools, BCFtools, GATK, STAR, salmon, MACS2) and still
power most genomics pipelines in 2026. But they also share a set of pain
points:

- **Memory unsafety.** Decades of CVEs, segfaults on edge-case BAMs, and bugs
  that only surface under threading. The C ecosystem has no answer for this
  other than "be careful."
- **Build fragility.** `./configure && make` against system zlib/htslib/curl,
  one toolchain per institution. Reproducibility relies on Conda papering
  over the cracks.
- **Concurrency afterthoughts.** Most tools added pthreads later; almost none
  scale linearly past 16 cores. Modern hardware is wasted.
- **Format bloat.** Each tool re-implements BAM/VCF parsing in subtly
  incompatible ways. htslib is the de-facto standard but is itself a C
  library with all the above issues.
- **No package manager.** A user assembling a pipeline manages a dozen
  unrelated build systems.

## Why Rust specifically

Rust is not the only memory-safe systems language, but for this domain it
ticks the boxes:

- **Memory safety without GC.** Important when streaming 100GB BAMs.
- **Zero-cost abstractions.** You can write the parser in idiomatic Rust and
  the compiled code is competitive with hand-tuned C.
- **`cargo`.** A real package manager. Pipelines become `cargo add`able.
- **SIMD and `rayon`.** First-class portable SIMD (stable as of 2024–25) and
  data parallelism that matches or beats OpenMP for most workloads.
- **PyO3 + Arrow + Polars.** Interop with the Python data-science stack is
  better than from C++, and shares a memory representation (Arrow) with the
  ecosystems we are *not* trying to displace.
- **Established beachhead.** `noodles`, `rust-bio`, `rust-htslib`, `needletail`,
  `minimap2-rs`, `alevin-fry`, `nf-core` all show Rust is viable in this
  domain. We are not pioneering — we are filling out a partial map.

## Why Zig, C++, etc. are not the choice

- **Zig** is promising and we acknowledge a few Zig bioinformatics projects
  (`htsz`, `zig-bio` experiments) but the ecosystem is two orders of
  magnitude smaller and the language is pre-1.0.
- **Modern C++** can match Rust on performance but not on package management
  or memory safety guarantees, and the existing C++ bioinformatics tools
  (e.g. samtools, BCFtools, STAR) demonstrate the problem rather than
  solving it.
- **Hand-written assembly** has a role for inner kernels (Smith-Waterman,
  hashing) but should be invoked from a safe host. Several existing Rust
  crates already FFI into vectorized C/asm — `wfa2-lib`, `ksw2`, etc.

## What success looks like

A user runs `cargo install rsomics-bwa rsomics-samtools rsomics-bcftools`,
gets statically-linked binaries that match or beat the C originals on a
standard benchmark, and writes a pipeline in 30 lines of Rust that previously
required Snakemake plus six Conda envs.

We are not there. This repo is the plan to get there.
