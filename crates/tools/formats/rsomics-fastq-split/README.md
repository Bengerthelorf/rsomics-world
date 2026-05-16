# rsomics-fastq-split

Split a FASTQ into many files — by file count or by line count (SE and PE).

```bash
cargo install rsomics-fastq-split
```

## Scope

The **split-only** partition of fastp's surface (one operation = one crate).
It partitions a large FASTQ into N pieces for parallel downstream processing.

| Operation | Crate |
|---|---|
| Split into N files / by line count | **rsomics-fastq-split** ← here |
| Adapter / poly-G / poly-X / fixed-length trim | rsomics-fastq-trim |
| Per-read quality + length filter | rsomics-fastq-filter |
| Inline-UMI extract + stamp | rsomics-fastq-umi |
| Exact / near dedup | rsomics-fastq-dedup |

## Behaviour

Output files are named `<zero-padded-index>.<basename(--out1)>` in `--out1`'s
directory (e.g. `out.fq.gz` → `0001.out.fq.gz`, `0002.out.fq.gz`, …), the
`.gz` suffix preserved so each file is independently gzip-compressed. In PE
mode R1/R2 are split in lockstep at the same record boundaries.

- **`--split_by_lines L`** — each file holds `L` lines (`L` must be a multiple
  of 4; a FASTQ record is 4 lines), the last file the remainder. Deterministic;
  **byte-equal to fastp 0.20.1** and verified so by the compat suite.
- **`--split N`** — exactly enough files of `ceil(total/N)` reads each, from a
  real read count. fastp's `--split N` distributes by a *file-size estimate*
  (not an exact count), so it is intentionally **not** a byte-compat target;
  ours is deterministic exact-count (a deliberate, documented improvement) and
  is gated by fastp-independent golden tests.

## Usage

```bash
# 8 exact-count files
rsomics-fastq-split -i in.fq.gz -o out.fq.gz --split 8

# 1000 reads (4000 lines) per file — fastp-0.20.1-equal
rsomics-fastq-split -i in.fq.gz -o out.fq.gz --split_by_lines 4000

# PE, JSON report
rsomics-fastq-split -i r1.fq.gz -I r2.fq.gz -o o1.fq.gz -O o2.fq.gz \
    --split 16 --json | jq .result
```

## Origin

Clean, independent Rust port informed by reading permissively-licensed
upstream source (allowed and cited):

- fastp split — `src/options.h` `SplitOptions` and the split writer, read at
  the `v0.20.1` tag (fastp, MIT). Paper: Chen et al. 2018, *Bioinformatics*,
  doi:10.1093/bioinformatics/bty560. `--split_by_lines` byte-faithful to
  0.20.1; `--split N` deliberately exact-count rather than fastp's estimate.
- seqkit `split2` — Shen et al. 2016, *PLOS ONE*, doi:10.1371/journal.pone.0163962
  (seqkit, MIT), secondary behavioural reference.

FASTQ reading is via `rsomics-seqio` (decode-only producer + parallel parse;
ISA-L igzip gz backend on Linux, pure-Rust flate2 elsewhere).

License: MIT OR Apache-2.0.
Upstream credit: [fastp](https://github.com/OpenGene/fastp) (MIT),
[seqkit](https://github.com/shenwei356/seqkit) (MIT).
