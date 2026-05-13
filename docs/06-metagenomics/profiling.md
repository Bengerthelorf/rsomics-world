# Metagenomic profiling

> Quantitative taxonomic and functional profiling: abundance re-estimation,
> pathway and gene-family inference, and microbial-load normalization.

## Scope

Includes: Bayesian/EM re-estimation of read assignments (Bracken),
pathway-level functional profiling (HUMAnN3), per-gene functional annotation
of metagenomic ORFs (eggNOG-mapper), average-genome-size and copy-number
normalization (MicrobeCensus), and KEGG/EggNOG ontology mapping. Excludes:
the underlying read classifier (see [classification](classification.md)) and
the MAG-resolved per-genome view (see [assembly-mag](assembly-mag.md)).

## Design notes

- HUMAnN3 is glue: MetaPhlAn4 for species pre-screen, ChocoPhlAn pangenome
  for nucleotide mapping, DIAMOND/UniRef for translated search, then SQL-ish
  joins on KEGG/MetaCyc. Most of the Python orchestration is replaceable by
  Rust + `polars`, but the value of any rewrite is bounded by the DIAMOND
  inner loop.
- Bracken is a few hundred lines of Python doing Bayesian re-estimation on
  Kraken outputs. Trivially portable to Rust once `rsomics-kraken` exists;
  a fast pure-Rust Bracken is mostly a polars-table merge.
- eggNOG-mapper is dominated by DIAMOND/MMseqs2 search time + a join against
  the eggNOG hierarchy. The orthology-transfer logic is the interesting
  Rust target; the search side defers to whatever protein-search crate we
  end up shipping.
- MicrobeCensus is small, focused, and would slot into `rsomics-meta` as a
  utility crate. The reference impl uses Python + HMMER/BLAST; pure Rust
  with `nthash` + streaming is straightforward.
- Most of these tools are I/O- and DB-bound, not CPU-bound — Rust's gains
  here are about pipeline ergonomics, memory discipline, and reproducible
  builds (cargo > pip + Conda for everything but the DB itself).
- License watch: HUMAnN3 is MIT, MetaPhlAn is MIT, Bracken is **GPL-3**,
  eggNOG-mapper is **GPL-3**, MicrobeCensus is **GPL-3**. Any direct port
  inherits GPL; clean-room rewrites do not.

## TODO

- [ ] **`HUMAnN3`** — community-wide functional pathway profiling.
  - Reference impl: `Python` (wraps Bowtie2 + DIAMOND + MetaPhlAn4) · [biobakery/humann](https://github.com/biobakery/humann) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `nHUMAnN` (academic reimplementation)
  - Parallelism: Python multiprocessing + downstream pthreads
  - SIMD: inherited from DIAMOND
  - Quadrant: —
  - GPU-amenable: maybe — protein search SIMT-friendly (DIAMOND analogues)
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-humann`)
  - Consumes primitives: `polars`, future `rsomics-diamond` or `rsomics-mmseqs`, `rsomics-metaphlan`, future `rsomics-stats`
  - Notes: Rust port = `polars` for the join/aggregation layer + delegate to a Rust DIAMOND-equivalent (open need) and to `rsomics-metaphlan`. The pathway-coverage math itself is the smallest part of the codebase.

- [ ] **`Bracken`** — Bayesian re-estimation of abundance from Kraken reports.
  - Reference impl: `Python` + `C` · [jenniferlu717/Bracken](https://github.com/jenniferlu717/Bracken) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: trivial parallel over taxa
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — table walking, latency-bound
  - Upstream license: `GPL-3`
  - Priority: `P0`
  - Layer: `subcommand-of-rsomics-kraken` (`--mode bracken` flag on the same binary)
  - Consumes primitives: `rsomics-kraken`, `polars`, `statrs`
  - Notes: Trivial port. Bracken essentially walks the Kraken report tree and redistributes reads by k-mer distribution priors. Once `rsomics-kraken` emits the same report format, Bracken is ~500 lines of Rust. Clean-room rewrite to avoid GPL inheritance.

- [ ] **`MicrobeCensus`** — average genome size + 16S rRNA copy-number normalization.
  - Reference impl: `Python` · [snayfach/MicrobeCensus](https://github.com/snayfach/MicrobeCensus) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — HMM scoring over universal genes parallelises
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-humann` (normalisation step inside the profiler umbrella)
  - Consumes primitives: `noodles-fastq`, `rsomics-kmer`, future `rsomics-hmm`
  - Notes: Aligns reads to ~30 universal single-copy gene models and divides. Small, well-bounded; useful normalizer for any downstream comparative analysis. Easy clean-room rewrite.

- [ ] **`MetaPhlAn4`** (profiling face) — quantitative species/strain abundance.
  - Reference impl: `Python` · [biobakery/MetaPhlAn](https://github.com/biobakery/MetaPhlAn) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `mOTUs`
  - Parallelism: Python multiprocessing
  - SIMD: inherited from Bowtie2
  - Quadrant: —
  - GPU-amenable: maybe — marker alignment SW-like
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-metaphlan` (the binary covers both classification and profiling — cross-listed with [`classification.md`](classification.md))
  - Consumes primitives: same as classification.md entry
  - Notes: Same crate as the classification entry. Profiling and classification are two faces of the same binary; do not split into two crates. **Cross-reference only — canonical entry is in `classification.md`.**

- [ ] **`Kaiju` + `functional` mode** — protein-space taxonomic + functional readout.
  - Reference impl: `C++` · [bioinformatics-centre/kaiju](https://github.com/bioinformatics-centre/kaiju) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — BWT protein-space search
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-mmseqs` (cross-listed with [`classification.md`](classification.md))
  - Consumes primitives: see classification.md
  - Notes: Kaiju emits a per-read protein hit list that can be joined to KEGG/COG. The interesting Rust work is the join, not the alignment. **Cross-reference only — canonical entry in `classification.md`.**

- [ ] **`eggNOG-mapper`** v2 — orthology-based functional annotation.
  - Reference impl: `Python` (wraps DIAMOND or MMseqs2 + HMMER) · [eggnogdb/eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper) · `GPL-3` (code) + `CC BY-NC` (database)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: inherited from DIAMOND
  - Quadrant: —
  - GPU-amenable: maybe — protein search SIMT-friendly
  - Upstream license: `GPL-3` (code); database `CC BY-NC`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-eggnog`)
  - Consumes primitives: future `rsomics-diamond` or `rsomics-mmseqs`, `polars`, future `rsomics-hmm`
  - Notes: The orthology-transfer + KEGG/GO/COG cross-walk is genuinely valuable. Database licensing (non-commercial) is a constraint on any bundled distribution. Treat code rewrite and DB packaging as separate problems.

- [ ] **`KEGG/GhostKOALA` and KO mappers**.
  - Reference impl: web service (KEGG) + open scripts · KEGG itself is closed-DB; tools like `kofamscan` are `BSD` · [takaram/kofam_scan](https://github.com/takaram/kofam_scan)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream serial
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — HMMER inner loops
  - Upstream license: `BSD` (kofamscan); KEGG DB closed
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-eggnog` (KO assignment mode)
  - Consumes primitives: future `rsomics-hmm`, `polars`
  - Notes: KEGG-API access is restricted (subscription). Focus on local HMMER-based KO assignment (`kofamscan`) which is open and re-implementable once we have a Rust HMMER-equivalent (separate open need).
