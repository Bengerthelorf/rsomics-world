# Conventions

Rules for contributing to this repo and for the Rust crates it ships.

## Repository shape

This is a single Cargo workspace, **not** a multi-repo setup.

```
rsomics-world/
├── Cargo.toml              workspace manifest
├── crates/
│   ├── foundation/         Layer A — library-only crates depended on by 2+ tools
│   │   └── rsomics-*/
│   └── tools/              Layer B — each is one installable binary
│       └── rsomics-*/
└── docs/                   per-domain planning, TODO checklists
```

- One git repo, one Cargo workspace, **independent semver per crate**.
- B → A → external. A never depends on B. B never depends on B; share via A.
- A is library-only (no `[[bin]]`). A crate with a binary is by definition B.
- Promote an internal module to A only when ≥ 2 B crates need it.

## Doc format

Every module sub-doc (the topic files under `docs/0X-*/`) follows this skeleton:

```markdown
# <Topic title>

> One-sentence elevator description.

## Scope

What this topic covers; what it does *not* cover (point to the neighboring
topic).

## Design notes

Algorithmic considerations, where Rust helps, where it doesn't. Brief — this
is not a research paper. 3–6 bullets is typical.

## TODO

A flat checklist. Each entry uses the schema below.
```

### TODO entry schema

```markdown
- [ ] **<canonical tool name>** — <one-line purpose>
  - Reference impl: <language> · <repo URL or "—"> · <license>
  - Existing Rust: <crate name + URL, or "none">
  - Existing Rust kind: pure-port | FFI-wrapper | partial-port | none
  - Existing non-C alternatives: <Zig / Go / C++ rewrites, or "—">
  - Parallelism: rayon | std::thread | tokio | async-runtime | single-threaded
  - SIMD: explicit (std::simd / packed_simd) | auto-vectorize | none
  - GPU-amenable: yes | no | maybe — <one-sentence why>
  - Upstream license: <SPDX>
  - Priority: P0 / P1 / P2
  - Layer: A (foundation) | B (tool) | adopt | subcommand-of-<crate>
  - Consumes primitives: <list of Layer A crates, or "—" for A entries>
  - Notes: <SIMD-critical? Already production-ready? Wrap vs rewrite? Algorithmic notes?>
```

Every tool gets an entry, *even if a mature Rust implementation already
exists*. Mark such entries with `[x]` and explain in the notes whether we
plan to (a) adopt as-is, (b) extend, or (c) leave alone.

### TODO.md legend

| Mark | Meaning |
|---|---|
| `[ ]` | Open — no Rust implementation, on our queue |
| `[~]` (FFI-wrapper) | Rust crate exists but wraps a C/C++ upstream via FFI |
| `[~]` (partial-port) | Pure-Rust port exists but is incomplete (stub, missing subcommand, single-threaded hot path) |
| `[x]` | Production-grade pure-Rust implementation exists — adopt |

## External-dependency quadrants

For every adopted external crate, classify and document its quadrant in the
TODO entry's `Existing Rust kind` and `Parallelism` / `SIMD` fields:

| Quadrant | What it is | Acceptable for | Note |
|---|---|---|---|
| ① Pure Rust + explicit parallelism (rayon / std::simd) | e.g. `noodles`, `rayon`, `polars`, `needletail` | Adopt freely | The goal state |
| ② FFI wrapper over C (`cc`, `bindgen`, `*-sys`) | e.g. `rust-htslib`, `minimap2-rs` | Adopt for speed; tag clearly as "FFI wrapper, not Rust rewrite" | C upstream stays the perf reference |
| ③ Pure Rust but single-threaded in hot path | varies — verify per case | **Avoid** in hot paths | Inherits the disease we're curing |
| ④ Pure Rust, non-hot, edge utility | `clap`, `serde`, `anyhow`, `regex` | Adopt freely | Doesn't matter |

Detection cues:
- Quadrant ②: `[build-dependencies] cc`/`bindgen`, `*-sys` crate name, presence of `build.rs`
- Quadrant ①: `rayon::`, `par_iter`, `std::simd::`, `packed_simd`, `target_feature`
- Quadrant ③: pure Rust with none of the above on the hot path
- Quadrant ④: small, edge-of-program crate

## License + clean-room

- This planning repo's docs: **CC BY 4.0**.
- All crates: dual **MIT OR Apache-2.0**.
- Crates derived from GPL upstream may stay MIT/Apache-2.0 only if a
  **clean-room methodology** is documented in that crate's README Origin
  section. The methodology must rely on the published paper, the public file
  format spec, and black-box behavior testing — never on reading the GPL
  source.

### Origin section template (GPL upstream)

Each crate that ports a GPL tool ships a README section like:

```markdown
## Origin

This crate is an independent Rust reimplementation of `<upstream>` based on:
- The published method (cite paper + DOI)
- The public file-format spec
- Black-box behavior testing against the upstream binary

No source code from the GPL upstream was read during implementation. Test
fixtures are independently generated or sourced from public benchmark
datasets (HG002 / 1000 Genomes / etc.).

License: MIT OR Apache-2.0.
Upstream credit: <upstream> <link> (<their license>).
```

Crates that wrap GPL upstream via FFI inherit the GPL link-time obligations;
those crates are flagged in the same Origin section with the explicit
license note.

## Cross-platform target

First-class — every CI commit runs on these four:

- `x86_64-unknown-linux-gnu`
- `aarch64-unknown-linux-gnu`
- `x86_64-apple-darwin`
- `aarch64-apple-darwin`

**Not supported**: Windows, other operating systems. Tools whose algorithm
hinges on Linux-only syscalls (e.g. `splice`) gate that path behind
`#[cfg(target_os = "linux")]` with a portable fallback.

## Test data tiers

Bench and integration tests draw from four tiers, sized for the environment
they run in:

| Tier | Location | Size budget | Use |
|---|---|---|---|
| 1 | inside git (`tests/golden/`) | < 100 MB total | Small synthetic FASTQ/BAM/VCF for unit tests |
| 2 | GHA runner (downloaded) | < 5 GB ephemeral | chr20/22 HG002 subsets, 1000 Genomes subsets — integration tests |
| 3 | local HDD (`/Volumes/Zane's HDD/rsomics-fixtures/`) | 50–200 GB | One real WGS sample, RNA-seq runs, smoke benchmarks |
| 4 | `ssh 4090:/data3/rsomics-fixtures/` | up to ~1 TB | Multi-sample cohort benches, large-scale perf validation |

Tier 1 lives in git. Tier 2 is downloaded on demand per
[`tests/fixtures-manifest.toml`](tests/fixtures-manifest.toml) (URL + SHA-256
+ size + license per fixture). Tier 3 and 4 are provisioned manually; tests
that need them read `BCMR_BENCH_DATA` and skip if unset.

## Crate naming

- Foundation crates: `rsomics-<primitive>` (e.g. `rsomics-common`,
  `rsomics-intervals`, `rsomics-kmer`, `rsomics-fm-index`, `rsomics-align-core`,
  `rsomics-stats`).
- Tool crates: `rsomics-<tool>` (e.g. `rsomics-fastp`, `rsomics-bam`,
  `rsomics-bwa`).
- A tool with multiple modes ships them as subcommands of a single binary
  (e.g. `rsomics-bam view`, `rsomics-bam sort`, `rsomics-bam markdup`), not as
  separate crates.

## Linking between docs

Use relative paths: `[topic](../02-genomics/alignment-short-read.md)`.
Cross-module references are encouraged — many tools straddle categories.

## Adding a new topic

1. Pick the module it belongs in. If it does not fit, raise an issue rather
   than inventing a 10th top-level module.
2. Add the file using the doc-format skeleton above.
3. Update the module's `README.md` index.
4. Add the entries to `TODO.md`.
