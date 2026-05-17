# rsomics-bam-flagstat

SAM/BAM flag statistics — count reads by SAM flag category with
QC-pass/fail split. Rust reimplementation of `samtools flagstat`.

## Install

```
cargo install rsomics-bam-flagstat
```

Single binary. Reads BAM (binary, indexed or unindexed), SAM (text),
or piped SAM from stdin (`-`).

## Usage

```
rsomics-bam-flagstat aligned.bam
rsomics-bam-flagstat -O json aligned.bam
samtools view -h in.cram | rsomics-bam-flagstat -
```

Output format matches `samtools flagstat` line-for-line (15 lines,
`QC-passed + QC-failed description (percentage)`) so downstream
parsers work unchanged. JSON mode (`-O json`) emits structured
counts for programmatic consumption.

## Origin

Independent Rust reimplementation of `samtools flagstat` based on:

- The SAM specification v1.6 (flag bit definitions, semantics).
- The published method: Li H et al., *The Sequence Alignment/Map
  format and SAMtools*, Bioinformatics 25:2078 (2009),
  [doi:10.1093/bioinformatics/btp352].
- Black-box behaviour testing against `samtools flagstat` (v1.19+).

samtools is MIT-licensed; the source was consulted for output format
and edge-case semantics (the diff-chr + mapQ>=5 line, primary vs
non-primary counting). BAM parsing uses [noodles] (pure Rust,
Quadrant ①).

License: MIT OR Apache-2.0.
Upstream credit: [samtools](https://github.com/samtools/samtools) (MIT).

### External-dep quadrant classification

- `noodles` (noodles-bam, noodles-sam) — Quadrant ① (pure Rust,
  spec-tracking SAM/BAM/CRAM library).
- `rsomics-common`, `rsomics-help`, `clap`, `serde`, `serde_json` —
  Quadrant ④ (edge utilities).

No FFI wrappers (no Quadrant ②); no single-threaded-in-hot-path deps
(no Quadrant ③). The hot path is a linear scan of BAM records —
noodles provides zero-copy lazy record decoding, so only the flag
u16 and a few optional fields are touched per record.

## Performance

The release contract: strictly faster wall-clock than `samtools
flagstat` on the perfgate fixture. `samtools flagstat` is already
fast (C + htslib's optimised BAM decoder); the Rust win comes from
noodles' lazy field decoding (flagstat only needs the flag u16 +
tid/mtid/mapq, not the full record) and zero-allocation iteration.

[noodles]: https://github.com/zaeleus/noodles
[doi:10.1093/bioinformatics/btp352]: https://doi.org/10.1093/bioinformatics/btp352
