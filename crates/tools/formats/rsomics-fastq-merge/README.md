# rsomics-fastq-merge

Merge overlapping paired-end FASTQ reads into single consensus reads —
a Rust port of `fastp`'s overlap merge (`fastp --merge`), with the
optional quality-aware overlap base correction (`fastp --correction`).

## Install

```
cargo install rsomics-fastq-merge
```

Single binary. R1/R2 gzip/bzip2/xz/zstd auto-detected via [needletail].

## Usage

```
rsomics-fastq-merge --in1 R1.fq.gz --in2 R2.fq.gz --correction -m merged.fq
rsomics-fastq-merge -i R1.fq -I R2.fq | gzip > merged.fq.gz
```

A pair is "overlapping" when read1 and the reverse complement of read2
share an alignment of at least `--overlap_len_require` (default 30) bp
with at most `--overlap_diff_limit` (5) mismatches and at most
`--overlap_diff_percent_limit` (20) % of the overlap mismatched —
fastp's exact two-pass detection (forward, then reverse for short
inserts / adapter read-through). The merged read is
`read1[..len1] + revcomp(read2)[overlap_len..]`, named
`<read1-name> merged_<len1>_<len2>`. With `--correction`, an in-overlap
mismatch where one mate is ≥Q30 and the other ≤Q14 is rewritten to the
high-quality mate's base before merging.

## Origin

Independent Rust reimplementation of the paired-end merge in
[`fastp`](https://github.com/OpenGene/fastp) (Chen et al., *fastp: an
ultra-fast all-in-one FASTQ preprocessor*, Bioinformatics 34:i884,
[doi:10.1093/bioinformatics/bty560]). fastp is MIT-licensed; the
algorithms and the default thresholds were read from the fastp **v0.20.1**
source (`src/overlapanalysis.cpp`, `src/basecorrector.cpp`,
`src/options.cpp`) and are cited rather than reconstructed. fastp's own
`overlapanalysis.cpp` / `basecorrector.cpp` `test()` vectors are
reproduced verbatim as unit-test oracles. The compat test is
version-gated to fastp v0.20.1 (newer fastp drifts; v0.20.1 is the
pinned reference installed in CI/publish).

License: MIT OR Apache-2.0. Upstream credit: [fastp] (MIT).

### External-dep quadrant classification

- `needletail` — Quadrant ① (pure Rust + SIMD).
- `rsomics-common`, `rsomics-help`, `clap`, `serde`, `serde_json`,
  `anyhow` — Quadrant ④ (edge utilities).

No FFI wrappers (no Quadrant ②); no known single-threaded-in-hot-path
deps (no Quadrant ③).

## Performance

The release contract: strictly faster wall-clock than `fastp --merge`
on the perfgate fixture (single thread, same input). Provenance lives
in `.autopilot/state/perf-*.md` and `benches/`.

[needletail]: https://crates.io/crates/needletail
[fastp]: https://github.com/OpenGene/fastp
[doi:10.1093/bioinformatics/bty560]: https://doi.org/10.1093/bioinformatics/bty560
