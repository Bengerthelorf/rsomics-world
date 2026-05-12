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
  with `nthash` + `rust-htslib`-style streaming is straightforward.
- Most of these tools are I/O- and DB-bound, not CPU-bound — Rust's gains
  here are about pipeline ergonomics, memory discipline, and reproducible
  builds (cargo > pip + Conda for everything but the DB itself).
- License watch: HUMAnN3 is MIT, MetaPhlAn is MIT, Bracken is **GPL-3**,
  eggNOG-mapper is **GPL-3**, MicrobeCensus is **GPL-3**. Any direct port
  inherits GPL; clean-room rewrites do not.

## TODO

- [ ] **`HUMAnN3`** — community-wide functional pathway profiling.
  - Reference impl: `Python` (wraps Bowtie2 + DIAMOND + MetaPhlAn4) · [biobakery/humann](https://github.com/biobakery/humann) · `MIT` (per setup.py)
  - Existing Rust: none verified
  - Existing non-C alternatives: `nHUMAnN` (academic reimplementation, check repo)
  - Priority: `P1`
  - Notes: Rust port = `polars` for the join/aggregation layer + delegate
    to a Rust DIAMOND-equivalent (open need) and to `rsomics-metaphlan`.
    The pathway-coverage math itself is the smallest part of the codebase.

- [ ] **`Bracken`** — Bayesian re-estimation of abundance from Kraken reports.
  - Reference impl: `Python` + `C` · [jenniferlu717/Bracken](https://github.com/jenniferlu717/Bracken) · `GPL-3`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P0` (bundled with Kraken2 work)
  - Notes: Trivial port. Bracken essentially walks the Kraken report tree
    and redistributes reads by k-mer distribution priors. Once
    `rsomics-kraken` emits the same report format, Bracken is ~500 lines
    of Rust. Clean-room rewrite to avoid GPL inheritance.

- [ ] **`MicrobeCensus`** — average genome size + 16S rRNA copy-number normalization.
  - Reference impl: `Python` · [snayfach/MicrobeCensus](https://github.com/snayfach/MicrobeCensus) · `GPL-3`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Aligns reads to ~30 universal single-copy gene models and divides.
    Small, well-bounded; useful normalizer for any downstream comparative
    analysis. Easy clean-room rewrite.

- [ ] **`MetaPhlAn4`** (profiling face) — quantitative species/strain abundance.
  - Reference impl: `Python` · [biobakery/MetaPhlAn](https://github.com/biobakery/MetaPhlAn) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: `mOTUs`
  - Priority: `P1`
  - Notes: Same crate as the classification entry. Profiling and classification
    are two faces of the same binary; do not split into two crates.

- [ ] **`Kaiju` + functional` mode** — protein-space taxonomic + functional readout.
  - Reference impl: `C++` · [bioinformatics-centre/kaiju](https://github.com/bioinformatics-centre/kaiju) · `GPL-3`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Kaiju emits a per-read protein hit list that can be joined to
    KEGG/COG. The interesting Rust work is the join, not the alignment.
    Shares an entry with the classification page.

- [ ] **`eggNOG-mapper`** v2 — orthology-based functional annotation.
  - Reference impl: `Python` (wraps DIAMOND or MMseqs2 + HMMER) · [eggnogdb/eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper) · `GPL-3` (code) + `CC BY-NC` (database)
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: The orthology-transfer + KEGG/GO/COG cross-walk is genuinely
    valuable. Database licensing (non-commercial) is a constraint on any
    bundled distribution. Treat code rewrite and DB packaging as separate
    problems.

- [ ] **`KEGG/GhostKOALA` and KO mappers**.
  - Reference impl: web service (KEGG) + open scripts · KEGG itself is closed-DB; tools like `kofamscan` are `BSD` · [taklam/KofamScan](https://github.com/takaram/kofam_scan)
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: KEGG-API access is restricted (subscription). Focus on local
    HMMER-based KO assignment (`kofamscan`) which is open and re-implementable
    once we have a Rust HMMER-equivalent (separate open need).
