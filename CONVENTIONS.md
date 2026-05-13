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
  - Existing Rust kind: pure-port | FFI-wrapper | partial-port | rust-native | none ( · slash-separated if multi-backend, dominant first, e.g. `pure-port/FFI-wrapper`)
  - Existing non-C alternatives: <Zig / Go / C++ rewrites, or "—">
  - Parallelism: rayon | std::thread | tokio | async-runtime | single-threaded
  - SIMD: explicit (std::simd / packed_simd) | auto-vectorize | none
  - Quadrant: ① | ② | ③ | ④
  - GPU-amenable: yes | maybe | no — <one-sentence judgment basis>
  - Upstream license: <SPDX>
  - Priority: P0 / P1 / P2
  - Layer: A (foundation) | B (tool) | adopt | subcommand-of-<crate>
  - Consumes primitives: <list of Layer A crates, or "—" for A entries>
  - Notes: <supplementary commentary; never the canonical answer for fields above>
```

Every tool gets an entry, *even if a mature Rust implementation already
exists*. Mark such entries with `[x]` and explain in the notes whether we
plan to (a) adopt as-is, (b) extend, or (c) leave alone.

#### `Existing Rust kind` — five values

- **`pure-port`** — Rust faithfully reimplements a specific C/C++ upstream (noodles ↔ htslib, needletail ↔ readfq, divsufsort ↔ libdivsufsort, lz4_flex ↔ lz4 C). Cross-validation against the upstream binary is the acceptance test.
- **`FFI-wrapper`** — Rust API over an existing C/C++ library (rust-htslib, libdeflater, hts-sys, cudarc, hdf5-metno-sys). Upstream C library stays the perf reference; license obligations follow the upstream's linking model on a per-crate basis.
- **`partial-port`** — Pure-Rust port covers a subset of the upstream's behavior (ruzstd decoder vs. zstd's full codec; rust-bio FM-index without FMD variant; AnnData crates missing `.uns`). Also covers **rust-native research-grade** Rust crates whose scope is partial relative to the canonical tool (the sub-case is disambiguated in `Notes`, not in the field value).
- **`rust-native`** — Rust-native concept with no specific upstream to port. Reference-impl field may cite an ecosystem analogue ("PyTorch", "Go runtime", "Intel TBB", "Mash") for orientation, but the Rust crate is not a code port. Examples: `rayon`, `tokio`, `candle`, `niffler`, `std::simd`, `fastbloom`, MinHash and HyperLogLog implementations of academic algorithms.
- **`none`** — no Rust implementation exists yet.

**Multi-backend convention.** When a single tool entry spans multiple backends with different kinds (e.g. flate2 with the default `miniz_oxide` is `pure-port`; with the `zlib-ng` feature is `FFI-wrapper`), write the field as a slash-separated list with the **dominant / default** backend first:

```
Existing Rust kind: pure-port/FFI-wrapper
```

The `Notes` line must state which backend is the default and which is feature-gated, so the reader can map slash-positions to backends.

#### `Quadrant` — canonical field, four values

- **①** — Pure-Rust with explicit parallelism (`rayon` / `std::simd` / `target_feature`) or pure-Rust orchestration over a perf-tuned codec. The goal state.
- **②** — FFI wrapper over C/C++ (`cc`, `bindgen`, `*-sys`, hand-FFI). C upstream stays the perf reference; adopt for speed but tag clearly as not a Rust rewrite.
- **③** — Pure-Rust but single-threaded in the hot path. **Avoid** in hot paths — inherits the disease we're curing. Surfaces opportunities for a Layer A replacement.
- **④** — Pure-Rust, non-hot, edge utility. Adoption decision doesn't move the needle on perf; pick on API/ergonomics. Examples: `clap`, `serde`, `anyhow`, `regex`, `niffler`-style format sniffers, small-format parsers.

A single entry may compose multiple quadrants across backends (e.g. `flate2` is ① with `miniz_oxide`, ② with `zlib-ng`). When that happens, write `Quadrant: ①+②` and let the Notes explain the split.

#### `GPU-amenable` — three values, strict definitions

- **`yes`** — clear GPU win; algorithm maps directly to SIMT (DL inference, large embarrassingly-parallel kernels, dense linear algebra). Speedup expected on relevant input sizes.
- **`maybe`** — requires algorithmic restructure to exploit GPU, OR only a portion of the work is GPU-amenable, OR the win is only at very large inputs. Engineering cost is non-trivial.
- **`no`** — inherently serial, I/O-bound, latency-dominated, or trivially small such that GPU offers no win (parsing, indexing, small lookups).

The trailing sentence on the line must state the judgment basis in one sentence ("dense matmul"; "irregular graph traversal"; "I/O-bound, parsing-only").

### TODO.md legend — four values

| Mark | Meaning |
|---|---|
| `[ ]` | Open — no Rust implementation, on our queue |
| `[~]` | A Rust crate exists, but is **incomplete in some dimension**: FFI-wrapper (wraps C/C++ via FFI), partial-port (covers a subset of upstream), OR rust-native but research-grade / partial-scope. The sub-case is named in the `Existing Rust kind` field and explained in `Notes`. |
| `[x]` | Production-grade pure-Rust implementation exists — **adopt the Rust crate** as a direct dependency. |
| `[A]` | **Adopt via subprocess** — the upstream is fine as-is (GPL-3 forbids linking; or hand-tuned SIMD that's not worth porting; or a stable Java binary that ships everywhere). We invoke it as a process from our pipelines but don't link or rewrite. Examples: Foldseek (GPL-3, hand-AVX), UShER (MIT but actively-maintained C++ at the scale we'd never match cheaply). |

The four marks are **mutually exclusive**. The `[A]` marker pairs with `Existing Rust kind: none` (we deliberately don't have a Rust crate; we adopt the upstream binary).

## External-dependency quadrants — detection cues

The `Quadrant` field is canonical (see the schema definitions above). This table gives the detection cues that drive the field value:

| Quadrant | Detection cue | Examples | Stance |
|---|---|---|---|
| ① pure-rust + explicit parallelism | `rayon::par_iter`, `std::simd`, `target_feature`, `crossbeam` channels for parallel pipelines | `noodles`, `rayon`, `polars`, `needletail`, `fastbloom`, `zarrs` | The goal state — adopt freely |
| ② FFI wrapper over C | `[build-dependencies] cc` or `bindgen`, `*-sys` crate, `build.rs` compiling C | `rust-htslib`, `libdeflater`, `cudarc`, `xz2`, `anndata-hdf5` | Adopt for speed; tag clearly. License obligations follow upstream's linking model per crate, not blanket |
| ③ pure-rust but single-threaded hot path | pure-Rust crate, no rayon / SIMD in its perf-critical modules | dated `rust-bio` FM-index path; some single-threaded parsers | Avoid in hot paths — opportunity for a Layer A replacement |
| ④ pure-rust edge utility | small, non-hot, edge-of-program | `clap`, `serde`, `anyhow`, `niffler` (orchestration), small parsers | Decision is API/ergonomics, not perf |

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

FFI wrappers inherit whatever license obligations the upstream's linking
model imposes — that is `htslib` (MIT/BSD-3) for `rust-htslib`, GPL-3.0 for
a hypothetical `bowtie2-sys`, etc. Document the case-by-case license note
in the wrapping crate's `## Origin` section. Never apply a blanket
inheritance assumption.

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
