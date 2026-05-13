# Mass spectrometry proteomics

> DDA and DIA peptide/protein identification, quantification, and FDR
> control from LC-MS/MS data.

## Scope

Includes: peptide-spectrum-match (PSM) search engines (Comet, MSFragger,
MSGF+, X!Tandem), integrated quantitative pipelines (MaxQuant, FragPipe,
DIA-NN), C++/Python/Java analysis frameworks (OpenMS), open data viewers
and quant tools (Skyline), and DIA-specific match-between-runs
(Quandenser, Spectronaut). Excludes: raw vendor format conversion (most
vendors ship closed-source converters; we'll wrap ProteoWizard
`msconvert` as FFI/process for the foreseeable future).

## Design notes

- The closed/proprietary footprint here is unusually large. MaxQuant
  (C#, free but not OSI), MSFragger (academic-only), DIA-NN (commercial
  from 1.9.2), Spectronaut (commercial). A pure-Rust open ecosystem
  exists in OpenMS (BSD-3, C++) and Comet (Apache-2.0); these are our
  natural friends.
- Core IO: mzML, mzXML, MGF, RAW (vendor). Rust has `mzdata` (and
  `mzpeaks`) which already cover mzML/MGF and IndexedMzML; build on it
  instead of reinventing. Vendor RAW conversion stays via `msconvert`.
- The scoring problem itself (PSM scoring against a target/decoy database)
  is well-bounded. Comet's ~10-15k lines of C++ are a tractable Rust port.
  MSFragger's fragment-ion indexing (the source of its speedup) is more
  involved but documented.
- DIA tools have a meaningful ML component (DIA-NN's neural-net peptide-
  property model, MSFragger's library predictions). Inference path is
  small models — `candle`/`burn` realistic.
- License watch: Comet **Apache-2.0**, X!Tandem **Artistic License 1.0**,
  OpenMS **BSD-3**, MSGF+ **proprietary-open**, MaxQuant **restricted-free**,
  MSFragger **academic-only**, FragPipe **academic-only**, DIA-NN **mixed/commercial**,
  Skyline **Apache-2.0**, Spectronaut **commercial**, Quandenser MIT/Apache.

## TODO

- [ ] **`MaxQuant`** — integrated DDA pipeline (search via Andromeda + MaxLFQ quant).
  - Reference impl: `C#` (.NET) · [maxquant.org](https://maxquant.org/) · proprietary-free (no public source)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `FragPipe`, `OpenMS`
  - Parallelism: .NET threading
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: no — closed source
  - Upstream license: proprietary-free
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-maxquant` as open MaxLFQ reimplementation on top of `Sage`)
  - Consumes primitives: `sage`, `mzdata`, `polars`, future `rsomics-stats`
  - Notes: Source-closed; only binaries available. No legal path to port. Rust strategy: implement open MaxLFQ-equivalent quant on top of an open search engine. The published MaxLFQ algorithm is reproducible.

- [ ] **`MSFragger`** — fragment-index-based ultra-fast PSM search.
  - Reference impl: `Java` · [Nesvilab/MSFragger](https://github.com/Nesvilab/MSFragger) · academic-only (not OSI)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `Comet`, `MSGF+`, `Sage` (Rust, see below)
  - Parallelism: JVM threading
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — fragment-ion matching is SIMT-trivial
  - Upstream license: academic-only
  - Priority: `P1`
  - Layer: —
  - Consumes primitives: —
  - Notes: Closed source, free for academic use only. Rust answer: **`Sage`** ([lazear/sage](https://github.com/lazear/sage)). Adopt Sage as the MSFragger-equivalent instead of porting.

- [x] **`Sage`** — MSFragger-style open Rust search engine.
  - Reference impl: `Rust` · [lazear/sage](https://github.com/lazear/sage) · `MIT`
  - Existing Rust: [`sage`](https://github.com/lazear/sage) (binary tool, install from source — crates.io name `sage` is squatted by an unrelated "age wrapper")
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: —
  - Parallelism: rayon
  - SIMD: auto-vectorize on fragment matching
  - Quadrant: ①
  - GPU-amenable: maybe — fragment indexing SIMT-friendly
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Adopt as the default search engine crate. Pure-Rust, fast, FDR-controlled. crates.io name squatted; install from source (squat catalog updated under module 07).

- [ ] **`Comet`** — open-source SEQUEST-style PSM search engine.
  - Reference impl: `C++` · [UWPR/Comet](https://github.com/UWPR/Comet) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `Sage` (Rust), `MSGF+` (Java)
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — same as MSFragger
  - Upstream license: `Apache-2.0`
  - Priority: `P2`
  - Layer: `subcommand-of-sage` (legacy-scoring mode on top of Sage if needed)
  - Consumes primitives: `mzdata`, `sage`
  - Notes: Sage mostly covers this need. Reimplementing Comet is only useful if we need its specific scoring quirks for legacy comparability.

- [ ] **`X!Tandem`** — early-2000s PSM search engine.
  - Reference impl: `C++` · [thegpm.org](http://www.thegpm.org/TANDEM/) · `Artistic License 1.0` — **canonical GitHub mirror at `thegpm/xtandem-vanilla` returns 404; project hosted on thegpm.org rather than GitHub**
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `Comet`, `MSGF+`, `Sage`
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: no — legacy
  - Upstream license: `Artistic License 1.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Legacy. Same logic as Comet — skip unless required for reproduction of older publications. Original entry's GitHub URL was wrong (404); the project lives on thegpm.org.

- [ ] **`MSGF+`** — Java search engine with universal scoring.
  - Reference impl: `Java` · [MSGFPlus/MSGFPlus](https://github.com/MSGFPlus/MSGFPlus) · BSD-style
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `Sage`, `Comet`
  - Parallelism: JVM threading
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — same as MSFragger
  - Upstream license: BSD-style
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Same comment as Comet/X!Tandem. Skip.

- [ ] **`OpenMS`** — comprehensive open MS framework.
  - Reference impl: `C++` · [OpenMS/OpenMS](https://github.com/OpenMS/OpenMS) · `BSD-3-Clause`
  - Existing Rust: none verified. `mzdata` provides an independent Rust mzML IO library that overlaps with one slice of OpenMS but is not a partial port — it's its own design
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `pyOpenMS` bindings
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — per-primitive
  - Upstream license: `BSD-3-Clause`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-ms-*` family of crates, each covering one primitive: peak picking, feature finding, isobaric quant)
  - Consumes primitives: `mzdata`, `mzpeaks`, `ndarray`, future `rsomics-stats`
  - Notes: Huge surface area (~1000+ tools). Don't reimplement wholesale. Identify a few high-value primitives and ship them as `rsomics-ms-*` crates. Defer the rest to OpenMS itself via FFI. `mzdata` provides a Rust mzML IO foundation we'd build on top of — not a Rust port of OpenMS.

- [ ] **`FragPipe`** — Java-based GUI/pipeline wrapping MSFragger + family.
  - Reference impl: `Java` · [Nesvilab/FragPipe](https://github.com/Nesvilab/FragPipe) · academic-only (mixed)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: JVM threading
  - SIMD: inherited from wrapped tools
  - Quadrant: —
  - GPU-amenable: no — orchestration
  - Upstream license: academic-only (mixed)
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Closed/restricted upstream. Rust analogue would be a thin CLI around `Sage` + open quant/FDR crates; FragPipe-feature parity is not a goal, just a usable pipeline.

- [ ] **`Skyline`** — targeted/DIA quant + visualization (open).
  - Reference impl: `C#` · [ProteoWizard/SkylineRunner](https://github.com/ProteoWizard) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: .NET threading
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: no — Windows-bound GUI
  - Upstream license: `Apache-2.0`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-ms-skyline-interop` (interop only, no GUI rewrite)
  - Consumes primitives: `mzdata`
  - Notes: Heavily Windows-bound GUI; little Rust value here. Provide interop instead — `.sky.zip` round-trip readers from `mzdata`.

- [ ] **`DIA-NN`** — DIA neural-net-based identification + quant.
  - Reference impl: `C++` · [vdemichev/DiaNN](https://github.com/vdemichev/DiaNN) · proprietary from 1.9.2 (versions ≤1.9.1 free)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `Sage`'s DIA mode (early-stage), `EncyclopeDIA`
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE/AVX
  - Quadrant: —
  - GPU-amenable: yes — neural-net RT/peptide property model
  - Upstream license: mixed (≤1.9.1 free; ≥1.9.2 commercial)
  - Priority: `P1`
  - Layer: `subcommand-of-sage` (Sage DIA mode)
  - Consumes primitives: `sage`, `candle` or `burn`, `mzdata`, future `rsomics-stats`
  - Notes: Most popular DIA tool, now commercial. Rust strategy: extend `Sage` with DIA mode + a small `candle`-based RT/peptide-property model. Real research problem, multi-year.

- [ ] **`Spectronaut`** (commercial, Biognosys) — DIA library-free search and quant.
  - Reference impl: closed-source · biognosys.com · commercial
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `DIA-NN`, `Sage` DIA-mode
  - Parallelism: closed
  - SIMD: closed
  - Quadrant: —
  - GPU-amenable: unknown — closed
  - Upstream license: commercial
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: No port path; closed. Listed for completeness.

- [ ] **`Quandenser`** — match-between-runs feature consolidation across DIA/DDA.
  - Reference impl: `C++` · [statisticalbiotechnology/quandenser](https://github.com/statisticalbiotechnology/quandenser) · `Apache-2.0` — **last update 2024-11-29, ~18 months stale at the threshold**
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — feature clustering parallelises
  - Upstream license: `Apache-2.0`
  - Priority: `P2`
  - Layer: `subcommand-of-sage` (MBR mode if needed)
  - Consumes primitives: `sage`, `mzdata`, `linfa-clustering`
  - Notes: Niche but well-bounded. C++ codebase is small (~10k lines). Reasonable Rust target if MBR becomes a priority. Upstream borderline-stale (similar to `paf` in module 01).

- [x] **`mzdata`** — pure-Rust mzML/MGF parser (proteomics IO foundation).
  - Reference impl: `Rust` · [mobiusklein/mzdata](https://github.com/mobiusklein/mzdata) · `Apache-2.0`
  - Existing Rust: [`mzdata`](https://crates.io/crates/mzdata) `0.63.5`; companion [`mzpeaks`](https://crates.io/crates/mzpeaks) `1.0.9`
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: —
  - Parallelism: rayon-able by caller
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: no — XML parsing, I/O bound
  - Upstream license: `Apache-2.0`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: The "noodles for proteomics". Adopt as the IO foundation; any new Rust MS tool should build on it.
