# rsomics-bbduk

K-mer-based contaminant removal and adapter/quality trimming for FASTQ —
an independent clean-room Rust reimplementation of **BBDuk** (BBTools,
JGI). "Duk" = Decontamination Using Kmers.

## Install

```
cargo install rsomics-bbduk
```

Single binary. Input gzip/bzip2/xz/zstd auto-detected via [needletail].

## Usage

```
# Right-trim Illumina adapters (full k-mer + partial 3'-tip via mink)
rsomics-bbduk -i R.fq.gz -r adapters.fa --ktrim r --mink 11 --hdist 1 -o clean.fq

# Drop PhiX-contaminated read pairs, divert them to a file
rsomics-bbduk -i R1.fq -I R2.fq -r phix.fa -o o1.fq -O o2.fq --outm phix.fq

# Literal reference + quality trim
rsomics-bbduk -i R.fq --literal AGATCGGAAGAGC --ktrim r --qtrim rl --trimq 10 -o c.fq
```

A reference is supplied as FASTA (`--ref`, repeatable) and/or inline
(`--literal a,b,c`). Every reference k-mer of length `--k` (default 27),
its reverse complement (unless `--no-rcomp`), and all k-mers within
`--hdist` (default 0) Hamming distance are indexed exactly. By default
the **middle base of every k-mer is a wildcard** (BBDuk `maskmiddle=t`,
disable with `--no-maskmiddle`); BBDuk turns this off automatically when
`--mink>0`. Short tip k-mers use their own Hamming distance `--hdist2`
(default 0 = exact tips, independent of `--hdist` — this matches BBDuk's
separate `hdist2` knob, not an omission).

- **`--ktrim f`** (default) — *kfilter*: a read sharing ≥ `--minkmerhits`
  (default 1) reference k-mers — or ≥ `--minkmerfraction` of its k-mers —
  is removed. Paired reads are kept or dropped **together**. Removed
  reads go to `--outm` if given.
- **`--ktrim r|l`** — trim the read from the first reference-k-mer hit to
  the 3' (`r`) or 5' (`l`) end. With `--mink N`, a partial adapter as
  short as `N` bp at the read **tip** is also trimmed (longest tip match
  wins).
- **`--qtrim r|l|rl`** then trims while the end base Phred is below
  `--trimq` (default 6).
- Reads shorter than `--minlength` (default 10) after trimming are
  discarded.

All defaults (`k=27`, `mink=0`, `hdist=0`, `hdist2=0`, `rcomp=t`,
`maskmiddle=t`, `mkh=1`, `mkf=0`, `ktrim=f`, `qtrim=f`, `trimq=6`,
`minlength=10`) match BBDuk's documented defaults so the default
invocation is byte-comparable to `bbduk.sh`.

`--k` is capped at 31 (the 2-bit codec limit; BBDuk's practical max) and
`--hdist` at 3 (beyond which the variant index outgrows memory) — both
fail loud rather than silently degrade.

## Origin

This crate is an independent Rust reimplementation of `BBDuk` (part of
BBTools / BBMap, Brian Bushnell, US DOE Joint Genome Institute) based on:

- The published method (DUK / BBDuk: Bushnell B., *BBMap: A Fast,
  Accurate, Splice-Aware Aligner*, LBNL-7065E; BBDuk documented
  algorithm: exact k-mer matching of reads against a reference k-mer set
  with reverse-complement, Hamming-distance (`hdist`) and shorter-tip
  (`mink`) extensions).
- The public `bbduk.sh` usage/help text and parameter documentation
  (default values, flag semantics).
- Black-box behaviour testing against the `bbduk.sh` binary.

**No BBTools/BBDuk source code was read or used as reference during
implementation.** BBTools is distributed under a restrictive
(non-OSI, redistribution-limited) license; this clean-room
methodology — published method + public CLI documentation + black-box
testing only — keeps this crate free of that license. Test fixtures are
independently generated.

License: MIT OR Apache-2.0.
Upstream credit: [BBTools/BBDuk](https://sourceforge.net/projects/bbmap/)
(Brian Bushnell, JGI; BBTools license — free, redistribution-restricted).

### External-dep quadrant classification

- `needletail` — Quadrant ① (pure Rust + SIMD).
- `rsomics-kmer` — Layer-A foundation (pure Rust 2-bit k-mer codec).
- `rsomics-common`, `rsomics-help`, `clap`, `serde`, `serde_json`,
  `anyhow` — Quadrant ④ (edge utilities).

No FFI wrappers (no Quadrant ②); no known single-threaded-in-hot-path
deps (no Quadrant ③). The reference index is an exact `HashSet<u64>` of
2-bit-encoded k-mers, **not** a Bloom filter: Bloom false positives
would flag clean reads and break byte-compatibility with BBDuk's
`kfilter`. A Bloom prefilter + exact-confirm is a documented future
bench-driven optimisation, not a correctness shortcut.

## Performance

The release contract: strictly faster wall-clock than `bbduk.sh` on the
perfgate fixture (single thread, same input, same defaults). Provenance
lives in `.autopilot/state/perf-*.md` and `benches/`.

[needletail]: https://crates.io/crates/needletail
