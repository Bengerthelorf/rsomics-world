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
  - Reference impl: `C#` (.NET) · [maxquant.org](https://maxquant.org/) · proprietary-free license (no public source)
  - Existing Rust: none verified
  - Existing non-C alternatives: `FragPipe`, `OpenMS` open alternatives
  - Priority: `P1`
  - Notes: Source-closed; only binaries available. No legal path to port.
    Rust strategy: implement open MaxLFQ-equivalent quant on top of an
    open search engine. The published MaxLFQ algorithm is reproducible.

- [ ] **`MSFragger`** — fragment-index-based ultra-fast PSM search.
  - Reference impl: `Java` · [Nesvilab/MSFragger](https://github.com/Nesvilab/MSFragger) · academic-only (not OSI)
  - Existing Rust: none verified
  - Existing non-C alternatives: `Comet`, `MSGF+`, `Sage` (Rust, see below)
  - Priority: `P1`
  - Notes: Closed source, free for academic use only. Rust answer: `Sage`
    ([lazear/sage](https://github.com/lazear/sage) ★) — a pure-Rust
    fragment-ion-indexing search engine inspired by MSFragger, MIT
    license, production-ready. **Adopt Sage** as the MSFragger-equivalent.

- [x] **`Sage`** — MSFragger-style open Rust search engine.
  - Reference impl: `Rust` · [lazear/sage](https://github.com/lazear/sage) · `MIT`
  - Existing Rust: same as reference
  - Existing non-C alternatives: —
  - Priority: `P0` (adopt)
  - Notes: Not in the original list but it is the Rust answer to MSFragger.
    Pure-Rust, fast, FDR-controlled. Adopt as the default search engine
    crate.

- [ ] **`Comet`** — open-source SEQUEST-style PSM search engine.
  - Reference impl: `C++` · [UWPR/Comet](https://github.com/UWPR/Comet) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: `Sage` (Rust), `MSGF+` (Java)
  - Priority: `P2`
  - Notes: `Sage` mostly covers this need. Reimplementing Comet is only
    useful if we need its specific scoring quirks for legacy comparability.

- [ ] **`X!Tandem`** — early-2000s PSM search engine, still cited.
  - Reference impl: `C++` · [thegpm/X-Tandem](https://github.com/thegpm/xtandem-vanilla) · `Artistic License 1.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: `Comet`, `MSGF+`, `Sage`
  - Priority: `P2`
  - Notes: Legacy. Same logic as Comet — skip unless required for
    reproduction of older publications.

- [ ] **`MSGF+`** — Java search engine with universal scoring.
  - Reference impl: `Java` · [MSGFPlus/MSGFPlus](https://github.com/MSGFPlus/MSGFPlus) · check repo (BSD-style for code; data restricted)
  - Existing Rust: none verified
  - Existing non-C alternatives: `Sage`, `Comet`
  - Priority: `P2`
  - Notes: Same comment as Comet/X!Tandem. Skip.

- [ ] **`OpenMS`** — comprehensive open MS framework.
  - Reference impl: `C++` · [OpenMS/OpenMS](https://github.com/OpenMS/OpenMS) · `BSD-3-Clause`
  - Existing Rust: none verified, but `mzdata` covers a slice
  - Existing non-C alternatives: `pyOpenMS` bindings; partial Rust via `mzdata`
  - Priority: `P1`
  - Notes: Huge surface area (~1000+ tools). Don't reimplement wholesale.
    Identify a few high-value primitives (peak picking, feature finding,
    isobaric quantitation) and ship them as `rsomics-ms-*` crates. Defer
    the rest to OpenMS itself via FFI.

- [ ] **`FragPipe`** — Java-based GUI/pipeline wrapping MSFragger + family.
  - Reference impl: `Java` · [Nesvilab/FragPipe](https://github.com/Nesvilab/FragPipe) · academic-only (mixed)
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Closed/restricted upstream. Rust analogue would be a thin CLI
    around `Sage` + open quant/FDR crates; FragPipe-feature parity is not
    a goal, just a usable pipeline.

- [ ] **`Skyline`** — targeted/DIA quant + visualization (open).
  - Reference impl: `C#` · [ProteoWizard/SkylineRunner](https://github.com/ProteoWizard) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Heavily Windows-bound GUI; little Rust value here. Provide
    interop instead — `.sky.zip` round-trip readers from `mzdata`.

- [ ] **`DIA-NN`** — DIA neural-net-based identification + quant.
  - Reference impl: `C++` · [vdemichev/DiaNN](https://github.com/vdemichev/DiaNN) · proprietary from 1.9.2 (versions ≤1.9.1 free)
  - Existing Rust: none verified
  - Existing non-C alternatives: `Sage`'s DIA mode (early-stage), `EncyclopeDIA`
  - Priority: `P1`
  - Notes: Most popular DIA tool, now commercial. Rust strategy: extend
    `Sage` with DIA mode + a small `candle`-based RT/peptide-property model.
    Real research problem, multi-year.

- [ ] **`Spectronaut`** (commercial, Biognosys) — DIA library-free search and quant.
  - Reference impl: closed-source · biognosys.com · commercial
  - Existing Rust: none verified
  - Existing non-C alternatives: `DIA-NN`, `Sage` DIA-mode
  - Priority: `P2`
  - Notes: No port path; closed. Listed for completeness.

- [ ] **`Quandenser`** — match-between-runs feature consolidation across DIA/DDA.
  - Reference impl: `C++` · [statisticalbiotechnology/quandenser](https://github.com/statisticalbiotechnology/quandenser) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Niche but well-bounded. C++ codebase is small (~10k lines).
    Reasonable Rust target if MBR becomes a priority.

- [x] **`mzdata`** — pure-Rust mzML/MGF parser (proteomics IO foundation).
  - Reference impl: `Rust` · [mobiusklein/mzdata](https://github.com/mobiusklein/mzdata) · `Apache-2.0`
  - Existing Rust: same as reference
  - Existing non-C alternatives: —
  - Priority: `P0` (adopt)
  - Notes: Listed even though not in original prompt, because it is the
    "noodles for proteomics". Adopt as the IO foundation; any new Rust
    MS tool should build on it.
