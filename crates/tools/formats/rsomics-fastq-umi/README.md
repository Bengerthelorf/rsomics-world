# rsomics-fastq-umi

Inline-UMI extraction and read-name stamping for FASTQ inputs (SE and PE).

```bash
cargo install rsomics-fastq-umi
```

## Scope

This crate is the **UMI-extract** partition of fastp's surface (one
operation = one crate; not a Swiss-army wrapper). It moves the 5' UMI off
the read and into the read name so downstream dedup/consensus tools can group
by molecule.

| Operation | Crate |
|---|---|
| UMI extract + stamp (all fastp `--umi_loc`) | **rsomics-fastq-umi** ← here |
| Adapter / poly-G / poly-X / fixed-length trim | rsomics-fastq-trim |
| Per-read quality + length filter | rsomics-fastq-filter |
| Exact / near dedup | rsomics-fastq-dedup |

The full fastp 0.20.1 `--umi_loc` set is implemented: `read1` / `read2`
(5' of that read's sequence), `index1` / `index2` (the read-name trailing
index field), `per_index`, `per_read`. `read2` and `index2` require PE.

## Behaviour

- **`read1` / `read2` / `per_read`** — `umi = seq[..min(umi_len, read_len)]`;
  the source read's seq+qual are trimmed at the 5' by `min(umi_len + umi_skip,
  read_len - 1)` (fastp `Read::trimFront` clamps to `length()-1`, keeping ≥1
  base). `per_read` (PE) merges both reads' UMIs as `umi1_umi2` and trims both.
- **`index1` / `index2` / `per_index`** — the UMI is the read-name index field
  (`firstIndex` of R1 / `lastIndex` of R2, fastp `src/read.cpp`); no trimming.
  `per_index` merges as `firstIndex_lastIndex`. A missing index field is a
  pass-through (fastp's `if(!umi.empty())` guard), not an error.

The UMI is written into the read name as `<delimiter>[<prefix>_]<umi>`,
inserted before the first space if the name has one, else appended. In PE both
mate names are stamped with the same UMI so the pair stays consistent. A
zero-length sequence-UMI source read fails loud (fastp's `trimFront` throws
there — no defined output; we error rather than fabricate one).

```text
@read1 1:N:0          ACGTACGTTTTT...   --umi_len 4 -->   @read1:ACGT 1:N:0   ACGTTTTT...
```

## Usage

```bash
# SE: 8 bp inline UMI on R1
rsomics-fastq-umi -i in.fq.gz -o out.fq.gz --umi_len 8

# PE: UMI on R1, JSON report
rsomics-fastq-umi -i r1.fq.gz -I r2.fq.gz -o o1.fq.gz -O o2.fq.gz \
    --umi_len 8 --json | jq .result
```

## Origin

This crate is a clean, independent Rust port. Implementation was informed by
reading permissively-licensed upstream source (allowed and cited):

- fastp UMI processing — `src/umiprocessor.cpp` (`UmiProcessor::process`,
  `addUmiToName`), `src/read.cpp` (`firstIndex`/`lastIndex`/`trimFront`), and
  `src/options.h` `UMIOptions`, read at the `v0.20.1` tag (fastp, MIT). Paper:
  Chen et al. 2018, *Bioinformatics*, doi:10.1093/bioinformatics/bty560.
- umi_tools `extract` — Smith, Heger, Sudbery 2017, *Genome Research*,
  doi:10.1101/gr.209601.116 (UMI-tools, MIT), secondary behavioural reference.

FASTQ reading is via `rsomics-seqio` (decode-only producer + parallel parse;
ISA-L igzip gz backend on Linux, pure-Rust flate2 elsewhere).

License: MIT OR Apache-2.0.
Upstream credit: [fastp](https://github.com/OpenGene/fastp) (MIT),
[UMI-tools](https://github.com/CGATOxford/UMI-tools) (MIT).
