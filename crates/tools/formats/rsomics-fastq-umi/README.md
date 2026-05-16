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
| Inline 5' UMI extract + stamp (`--umi_loc read1`/`read2`) | **rsomics-fastq-umi** ← here |
| Adapter / poly-G / poly-X / fixed-length trim | rsomics-fastq-trim |
| Per-read quality + length filter | rsomics-fastq-filter |
| Exact / near dedup | rsomics-fastq-dedup |

Index-based UMI locations (`index1`/`index2`/`per_index`/`per_read`) are not
implemented — they need separate index-FASTQ inputs or dual-UMI merge and no
consumer needs them yet (YAGNI). The inline read1/read2 case is the dominant
UMI workflow and the `umi_tools extract` default.

## Behaviour

For the UMI-carrying read: `umi = seq[..min(umi_len, read_len)]`; the read's
sequence and quality are then trimmed at the 5' by `min(umi_len, read_len) +
umi_skip` (fastp `trimFront(umiLength + umi_skip)`). The UMI is written into
the read name as `<delimiter>[<prefix>_]<umi>`, inserted before the first
space if the name has one, else appended. In PE, only the UMI-source read is
trimmed; both mate names are stamped with the same UMI so the pair stays
consistent.

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
  `addUmiToName`) and `src/options.h` `UMIOptions` (fastp, MIT). Paper:
  Chen et al. 2018, *Bioinformatics*, doi:10.1093/bioinformatics/bty560.
- umi_tools `extract` — Smith, Heger, Sudbery 2017, *Genome Research*,
  doi:10.1101/gr.209601.116 (UMI-tools, MIT), secondary behavioural reference.

FASTQ reading is via `rsomics-seqio` (decode-only producer + parallel parse;
ISA-L igzip gz backend on Linux, pure-Rust flate2 elsewhere).

License: MIT OR Apache-2.0.
Upstream credit: [fastp](https://github.com/OpenGene/fastp) (MIT),
[UMI-tools](https://github.com/CGATOxford/UMI-tools) (MIT).
