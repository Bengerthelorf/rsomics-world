# Preprocessing & QC

> Quality trimming, adapter removal, FASTQ wrangling, and read-level QC
> before alignment.

## Scope

Everything that touches raw FASTQ before the aligner sees it: trimmers
(fastp, Trimmomatic, cutadapt, BBDuk, Trim Galore, AdapterRemoval), QC
reporters (FastQC, MultiQC), and Swiss-Army-knife sequence-manipulation
tools (seqkit, seqtk). Excludes basecalling (Guppy/Dorado, vendor-specific
and largely closed-source) and post-alignment QC (samtools stats, mosdepth
— covered under module 09 utility).

## Design notes

- `fastp` is the runaway winner in 2026 production pipelines: ultra-fast
  C/C++, all-in-one (trim + filter + QC + UMI), JSON+HTML reports, MultiQC-
  compatible. A pure-Rust replacement (`fasterp`) already exists and aims
  for bit-identical JSON output.
- The combined-tool advantage of fastp (one pass over the data) is the
  performance lever for any Rust rewrite. Avoid the
  cutadapt-then-trimmomatic-then-fastqc pattern; ship a single binary
  that does it all.
- For QC, RastQC and fastqc-rs are credible Rust FastQC alternatives.
  MultiQC remains the report-aggregation layer; a Rust port is *not* a
  priority (Python is fine for report rendering).
- `seqkit` (Go) and `seqtk` (C) are the universal sequence-manipulation
  CLIs; both have Rust ports (`seqtk-rs`, `faster`). A unified
  `rsomics-seqtools` umbrella matching the seqkit subcommand surface is
  worth the investment.
- Per-domain trimmers (AdapterRemoval for ancient DNA, BBDuk for kmer-
  based contaminant removal) have unique value; they are not redundant
  with fastp.

## TODO

- [x] **`fastp`** — all-in-one FASTQ preprocessor (Rust port exists).
  - Reference impl: `C/C++` · [OpenGene/fastp](https://github.com/OpenGene/fastp) · `MIT`
  - Existing Rust: [`fasterp`](https://docs.rs/fasterp) (pure-Rust port,
    fastp-compatible JSON output)
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Adopt fasterp as the target. Audit feature parity (UMI,
    polyG trimming, overlap correction, output splitting). One-pass
    SIMD-friendly inner loops on quality + base counting.

- [ ] **`Trimmomatic`** — adapter + quality trimming (Java).
  - Reference impl: `Java` · [usadellab/Trimmomatic](https://github.com/usadellab/Trimmomatic) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Largely superseded by fastp on speed and feature set. GPL
    licence. Document interop only — many legacy pipelines still call it.

- [ ] **`cutadapt`** — gold-standard adapter-trimming tool.
  - Reference impl: `Python / C (Cython)` · [marcelm/cutadapt](https://github.com/marcelm/cutadapt) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Reference implementation for adapter alignment semantics.
    `fastp` and `fasterp` implement adapter trimming differently
    (overlap-based) — for protocols that need cutadapt's specific
    semantics (small RNA, single-cell linker handling) a Rust port
    matters. MIT licence.

- [ ] **`BBDuk`** — k-mer-based contaminant removal + trimmer.
  - Reference impl: `Java` · [BBTools (JGI)](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) · proprietary-but-free
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Unique among trimmers for k-mer contaminant matching against
    a reference (phiX, rRNA, host). The Java is slow and memory-hungry;
    a Rust port using `fastbloom` + `nthash` for the k-mer side would
    be a clear win.

- [ ] **`Trim Galore`** — cutadapt + FastQC wrapper (Perl).
  - Reference impl: `Perl` · [FelixKrueger/TrimGalore](https://github.com/FelixKrueger/TrimGalore) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Effectively superseded by fastp. GPL-3.0. Document only.

- [ ] **`AdapterRemoval`** — paired-end + ancient-DNA trimmer.
  - Reference impl: `C++` · [MikkelSchubert/adapterremoval](https://github.com/MikkelSchubert/adapterremoval) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Distinctive value for ancient-DNA workflows (collapse paired
    reads, deal with short fragments). GPL-3.0. Lower priority than
    fastp.

- [~] **`FastQC`** — per-base QC report generator.
  - Reference impl: `Java` · [s-andrews/FastQC](https://github.com/s-andrews/FastQC) · `GPL-3.0`
  - Existing Rust: [`fastqc-rs`](https://fastqc-rs.github.io/);
    [`RastQC`](https://github.com/Huang-lab/RastQC) (2-3× faster,
    streaming, single static binary)
  - Existing non-C alternatives: `falco` (C++ port);
    `seqkit stats` (Go)
  - Priority: `P0`
  - Notes: RastQC is the strongest pure-Rust candidate; adopt and extend.
    MultiQC-compatible output is a hard requirement.

- [ ] **`MultiQC`** — report aggregator.
  - Reference impl: `Python` · [MultiQC/MultiQC](https://github.com/MultiQC/MultiQC) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Pure aggregation + templating. Python is *fine* for report
    rendering. Document the output schema and make sure every Rust
    tool (`fasterp`, `RastQC`, `varlociraptor`, etc.) emits MultiQC-
    parseable JSON. No Rust rewrite needed.

- [~] **`seqkit`** — Swiss-Army-knife FASTA/FASTQ tool (Go).
  - Reference impl: `Go` · [shenwei356/seqkit](https://github.com/shenwei356/seqkit) · `MIT`
  - Existing Rust: [`faster`](https://github.com/angelovangel/faster)
    (partial, fastq stats);
    [`seqtk-rs`](https://github.com/yenyen1/seqtk-rs) (seqtk-like subset)
  - Existing non-C alternatives: `seqkit` itself (Go); pyfastx (Python)
  - Priority: `P1`
  - Notes: Go performance is already strong; the win for Rust is
    deeper integration with `noodles` and the rest of the rsomics
    stack. A `rsomics-seqkit` matching the most-used 80% of seqkit
    subcommands is the right scope.

- [~] **`seqtk`** — lh3's minimalist sequence toolkit.
  - Reference impl: `C` · [lh3/seqtk](https://github.com/lh3/seqtk) · `MIT`
  - Existing Rust: [`seqtk-rs`](https://github.com/yenyen1/seqtk-rs)
  - Existing non-C alternatives: `seqkit` (Go, superset functionality)
  - Priority: `P1`
  - Notes: `seqtk-rs` covers the common subcommands. Adopt as a
    starting point for the broader `rsomics-seqkit` work.
