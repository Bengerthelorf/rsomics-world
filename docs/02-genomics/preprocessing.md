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

- [x] **`fastp`** — all-in-one FASTQ preprocessor.
  - Reference impl: `C/C++` · [OpenGene/fastp](https://github.com/OpenGene/fastp) · `MIT`
  - Existing Rust: [`fasterp`](https://crates.io/crates/fasterp) `0.2.1` (pure-Rust port, fastp-compatible JSON output; default `native` feature pulls `libdeflater` + `flate2/zlib-ng` FFI codec, plus AWS SDK + tokio for cloud streaming)
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: —
  - Parallelism: rayon over read pairs
  - SIMD: auto-vectorize on quality + base counting; explicit candidates for `std::simd` in adapter overlap scan; codec SIMD inherited from `libdeflater` / `zlib-ng`
  - Quadrant: ①+②
  - GPU-amenable: maybe — read-level parallelism is SIMT-trivial; engineering cost vs. CPU performance is marginal
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-fastp`, written original from scratch rather than adopting fasterp wholesale; see Notes)
  - Consumes primitives: `needletail`, `noodles-fastq`, future `rsomics-stats`
  - Notes: Audit feature parity (UMI, polyG trimming, overlap correction, output splitting). One-pass SIMD-friendly inner loops on quality + base counting. **Phase 2 verdict**: write `rsomics-fastp` original rather than adopt fasterp. Reasons: (1) fasterp's `default = ["native"]` feature pulls FFI codec backends (`libdeflater`, `flate2/zlib-ng`), making its hot path Quadrant ② not the ① this entry originally claimed; (2) fasterp's AWS-SDK + tokio dep surface is wider than a focused fastp-equivalent needs; (3) writing original lets us own the codec backend choice (zlib-rs default; libdeflater behind a feature flag) and establish the rsomics-* tool playbook (compat.rs + benches + golden tests) on a clean slate. fasterp + upstream fastp serve as compat-test references.

- [ ] **`Trimmomatic`** — adapter + quality trimming (Java).
  - Reference impl: `Java` · [usadellab/Trimmomatic](https://github.com/usadellab/Trimmomatic) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: JVM threading
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — same constraints as fastp
  - Upstream license: `GPL-3.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Largely superseded by fastp on speed and feature set. GPL licence. Document interop only — many legacy pipelines still call it.

- [ ] **`cutadapt`** — gold-standard adapter-trimming tool.
  - Reference impl: `Python / C (Cython)` · [marcelm/cutadapt](https://github.com/marcelm/cutadapt) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: upstream's Cython adapter alignment uses SSE
  - Quadrant: —
  - GPU-amenable: maybe — adapter alignment is SW; per-read SIMT-trivial
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-fastp` (cutadapt's specific adapter semantics surface as a `--cutadapt-compat` mode rather than a separate binary)
  - Consumes primitives: `needletail`, `noodles-fastq`, `block-aligner` for adapter SW
  - Notes: Reference implementation for adapter alignment semantics. `fastp` and `fasterp` implement adapter trimming differently (overlap-based) — for protocols that need cutadapt's specific semantics (small RNA, single-cell linker handling) a Rust port matters. MIT licence.

- [ ] **`BBDuk`** — k-mer-based contaminant removal + trimmer.
  - Reference impl: `Java` · [BBTools (JGI)](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) · proprietary-but-free
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: JVM threading
  - SIMD: none significant
  - Quadrant: —
  - GPU-amenable: maybe — k-mer hashing is SIMT-trivial
  - Upstream license: proprietary-but-free
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-bbduk`)
  - Consumes primitives: `rsomics-kmer` (`nthash-rs`), `fastbloom`, `needletail`, `noodles-fastq`
  - Notes: Unique among trimmers for k-mer contaminant matching against a reference (phiX, rRNA, host). The Java is slow and memory-hungry; a Rust port using `fastbloom` + `nthash` for the k-mer side would be a clear win.

- [ ] **`BFC`** — k-mer-spectrum substitution-error corrector for Illumina reads.
  - Reference impl: `C` · [lh3/bfc](https://github.com/lh3/bfc) · `MIT`
  - Existing Rust: none verified (`cargo search "error correction fastq"` → only filter/simulate/split-kmer crates; no read corrector)
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `Lighter` ([mourisl/Lighter](https://github.com/mourisl/Lighter), C++, `GPL-3.0` — Bloom-cardinality corrector; clean-room only); `Karect`, `Musket`, `Bloocoo` (older, lower priority)
  - Parallelism: BFC pthreads over reads; ours rayon over read batches
  - SIMD: k-mer roll (ntHash) auto-vectorises; correction-decision path is branchy
  - Quadrant: ① (pure Rust + rayon; `rsomics-kmer` ntHash + `fastbloom` blocked filter)
  - GPU-amenable: maybe — k-mer counting/lookup is SIMT-trivial; the greedy correction walk is branch-divergent, CPU is the right first target
  - Upstream license: `MIT` (read + cite `bfc.c` permitted; it is the behavioural + compat reference)
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-fastq-correct`; one operation = correct substitution errors by k-mer spectrum; a standalone pre-assembly / pre-variant-calling pipeline stage)
  - Consumes primitives: `rsomics-kmer` (k-mer count + ntHash, Layer-A 0.1.0), `rsomics-seqio` (FASTQ in), `rsomics-fqgz` (gz out), `fastbloom`
  - Notes: BFC builds a k-mer coverage histogram, picks a trusted/untrusted threshold from the bimodal spectrum, then greedily rewrites untrusted bases along the highest-support k-mer path. Write original (BFC is the MIT compat oracle; Lighter GPL → paper + black-box only). Feature-complete = BFC's real flag surface (`-k` k-mer, `-s` genome-size/auto, `-t` threads, drop-vs-correct policy), not an MVP. golden + version-gated compat vs `bfc`; 4090 perfgate `>1.0×`; L2 FreshEye on the threshold-pick + correction-walk fidelity to `bfc.c`.

- [ ] **`Trim Galore`** — cutadapt + FastQC wrapper (Perl).
  - Reference impl: `Perl` · [FelixKrueger/TrimGalore](https://github.com/FelixKrueger/TrimGalore) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Perl shell-level
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — orchestration wrapper
  - Upstream license: `GPL-3.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Effectively superseded by fastp. GPL-3.0. Document only.

- [ ] **`AdapterRemoval`** — paired-end + ancient-DNA trimmer.
  - Reference impl: `C++` · [MikkelSchubert/adapterremoval](https://github.com/MikkelSchubert/adapterremoval) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: no — same constraints as fastp
  - Upstream license: `GPL-3.0`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-fastp` (ancient-DNA / collapse mode)
  - Consumes primitives: `needletail`, `noodles-fastq`, `block-aligner`
  - Notes: Distinctive value for ancient-DNA workflows (collapse paired reads, deal with short fragments). GPL-3.0. Lower priority than fastp.

- [~] **`FastQC`** — per-base QC report generator.
  - Reference impl: `Java` · [s-andrews/FastQC](https://github.com/s-andrews/FastQC) · `GPL-3.0`
  - Existing Rust: [`fastqc-rs`](https://crates.io/crates/fastqc-rs) `0.3.4` ([fastqc-rs/fastqc-rs](https://github.com/fastqc-rs/fastqc-rs)); [`RastQC`](https://github.com/Huang-lab/RastQC) (binary tool, install from source, 2-3× faster, streaming, single static binary)
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: `falco` (C++ port); `seqkit stats` (Go)
  - Parallelism: fastqc-rs is rayon-able per-chunk; RastQC streaming
  - SIMD: auto-vectorize on quality histograms
  - Quadrant: ①
  - GPU-amenable: no — histogram accumulation is memory-bandwidth-bound
  - Upstream license: `GPL-3.0`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-fastqc` based on RastQC)
  - Consumes primitives: `needletail`, `noodles-fastq`, `rsomics-stats`
  - Notes: RastQC is the strongest pure-Rust candidate; adopt and extend. MultiQC-compatible output is a hard requirement.

- [ ] **`MultiQC`** — report aggregator.
  - Reference impl: `Python` · [MultiQC/MultiQC](https://github.com/MultiQC/MultiQC) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python serial
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — report aggregation
  - Upstream license: `GPL-3.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Pure aggregation + templating. Python is *fine* for report rendering. Document the output schema and make sure every Rust tool (`fasterp`, `RastQC`, `varlociraptor`, etc.) emits MultiQC-parseable JSON. No Rust rewrite needed.

- [~] **`seqkit`** — Swiss-Army-knife FASTA/FASTQ tool (Go).
  - Reference impl: `Go` · [shenwei356/seqkit](https://github.com/shenwei356/seqkit) · `MIT`
  - Existing Rust: [`faster`](https://github.com/angelovangel/faster) (binary tool, install from source — crates.io name `faster` is squatted by an unrelated explicit-SIMD library); [`seqtk-rs`](https://crates.io/crates/seqtk-rs) `0.2.0` (seqtk-like subset)
  - Existing Rust kind: `partial-port`
  - Existing non-C alternatives: `seqkit` itself (Go); `pyfastx` (Python)
  - Parallelism: rayon (faster); single-threaded (seqtk-rs)
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: no — text-based sequence manipulation
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-seqkit`)
  - Consumes primitives: `needletail`, `noodles-fasta`, `noodles-fastq`
  - Notes: Go performance is already strong; the win for Rust is deeper integration with `noodles` and the rest of the rsomics stack. A `rsomics-seqkit` matching the most-used 80% of seqkit subcommands is the right scope.

- [~] **`seqtk`** — lh3's minimalist sequence toolkit.
  - Reference impl: `C` · [lh3/seqtk](https://github.com/lh3/seqtk) · `MIT`
  - Existing Rust: [`seqtk-rs`](https://crates.io/crates/seqtk-rs) `0.2.0`
  - Existing Rust kind: `partial-port`
  - Existing non-C alternatives: `seqkit` (Go, superset functionality)
  - Parallelism: single-threaded today
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: no — text-based sequence manipulation
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-seqkit` (single binary; seqtk subset of subcommands)
  - Consumes primitives: `needletail`, `noodles-fasta`, `noodles-fastq`
  - Notes: `seqtk-rs` covers the common subcommands. Adopt as a starting point for the broader `rsomics-seqkit` work.
