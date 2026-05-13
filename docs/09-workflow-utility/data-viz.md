# Data visualization

> Genome browsers, plot libraries, tree viewers, and assembly-graph viewers.

## Scope

Includes: interactive genome browsers (IGV, JBrowse2, UCSC Browser, Gosling),
bulk-plot generators (deepTools, ggplot2, Plotly), tree visualization
(ETE3, iTOL), assembly-graph viewers (Bandage), and domain-specific
viewers (gnomAD viewer). Excludes: structure visualization (covered in
[../07-proteomics-structure/structure-analysis.md](../07-proteomics-structure/structure-analysis.md)).

## Design notes

- The visualization layer is dominated by **JavaScript** (IGV.js, JBrowse2,
  Gosling, igv-reports) and **R** (ggplot2). Rust's play is mostly back-end.
- For tree viz (ETE3, iTOL), Rust has no contender. Provide good Newick /
  NEXUS IO via a canonical Rust tree crate.
- Bandage is Qt/C++ — a desktop GUI. No Rust play.
- License watch: IGV **MIT**, JBrowse2 **Apache-2.0**, UCSC Genome
  Browser **non-commercial source-available**, Gosling **MIT**, deepTools
  **GPL-3**, ggplot2 **MIT**, Plotly.js **MIT**, ETE3 **GPL-3**, Bandage
  **GPL-3**.

## TODO

- [ ] **`IGV` / `IGV.js`** — Integrative Genomics Viewer (desktop + web).
  - Reference impl: `Java` (desktop) + `JavaScript` (web) · [igvteam/igv](https://github.com/igvteam/igv) and [igvteam/igv.js](https://github.com/igvteam/igv.js) · `MIT`
  - Existing Rust: none for viewer itself; backend file-serving via `bigtools`, `noodles`
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `JBrowse2`
  - Parallelism: JVM / browser
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — viewer
  - Upstream license: `MIT`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: `bigtools`, `noodles` (as data backends)
  - Notes: Don't rewrite the viewer. Provide solid Rust IO so IGV can load our outputs (CRAM, BigWig, VCF) over HTTP range requests.

- [ ] **`JBrowse2`** — modern JavaScript / TypeScript genome browser.
  - Reference impl: `TypeScript` · [GMOD/jbrowse-components](https://github.com/GMOD/jbrowse-components) · `Apache-2.0`
  - Existing Rust: none for viewer; WASM-based plugin path exists in JBrowse2
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: browser
  - SIMD: WebAssembly SIMD potentially
  - Quadrant: —
  - GPU-amenable: maybe — WebGL rendering
  - Upstream license: `Apache-2.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: future `rsomics-wasm` plugins
  - Notes: WASM-compiled Rust plugins for JBrowse2 are a realistic Phase 5+ direction. Until then, focus on backend.

- [ ] **`UCSC Genome Browser`** — classic web genome browser + Kent tools.
  - Reference impl: `C` · [ucscGenomeBrowser/kent](https://github.com/ucscGenomeBrowser/kent) · UCSC source-available (non-commercial)
  - Existing Rust: [`bigtools`](https://crates.io/crates/bigtools) `0.5.6` — pure-Rust BigWig/BigBed reader/writer (replaces parts of `kent`)
  - Existing Rust kind: `partial-port` (covers BigWig/BigBed slice of kent)
  - Existing non-C alternatives: —
  - Parallelism: rayon (bigtools)
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: no — file format IO
  - Upstream license: UCSC source-available
  - Priority: `P1`
  - Layer: `adopt` (bigtools)
  - Consumes primitives: —
  - Notes: We do not rebuild the browser. The Rust win is `bigtools` and related format readers/writers, which let pipelines avoid the UCSC license entirely. Adopt `bigtools`.

- [x] **`bigtools`** — Rust BigWig/BigBed library + CLI.
  - Reference impl: `Rust` · [jackh726/bigtools](https://github.com/jackh726/bigtools) · `MIT`
  - Existing Rust: [`bigtools`](https://crates.io/crates/bigtools) `0.5.6`
  - Existing Rust kind: `rust-native` (algorithm is the crate's own contribution; BigWig/BigBed spec is independent from kent code)
  - Existing non-C alternatives: —
  - Parallelism: rayon
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: no — file format IO
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Foundational. Use it wherever pipelines need BigWig/BigBed.

- [ ] **`Gosling`** — grammar-based interactive genomic viz.
  - Reference impl: `TypeScript` · [gosling-lang/gosling.js](https://github.com/gosling-lang/gosling.js) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: browser
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — browser viz
  - Upstream license: `MIT`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: JS only; Rust play is producing Gosling-spec-compatible JSON output from analysis pipelines.

- [ ] **`deepTools`** — Python coverage / heatmap plotting from BAM/BigWig.
  - Reference impl: `Python` · [deeptools/deepTools](https://github.com/deeptools/deepTools) · `GPL-3`
  - Existing Rust: none verified (but the pieces exist — `noodles-bam`, `bigtools`, `plotters`)
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: maybe — coverage aggregation parallelises trivially
  - Upstream license: `GPL-3`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-deeptools`)
  - Consumes primitives: `noodles-bam`, `bigtools`, `rsomics-coverage`, `plotters` or `plotly`, `rsomics-intervals`, `rayon`
  - Notes: Excellent Rust rewrite target. Inner loops (coverage aggregation over millions of intervals × many BAMs) are CPU-bound and parallelize trivially. A `rsomics-deeptools`-equivalent would be measurably faster and avoid Python deps for pipelines.

- [ ] **`ggplot2`** — R grammar-of-graphics plotting (de facto standard).
  - Reference impl: `R` · [tidyverse/ggplot2](https://github.com/tidyverse/ggplot2) · `MIT`
  - Existing Rust: partial — [`plotters`](https://crates.io/crates/plotters) `0.3.7` ([`plotters-rs/plotters`](https://github.com/plotters-rs/plotters)), [`plotly`](https://crates.io/crates/plotly) `0.14.1` ([`plotly/plotly.rs`](https://github.com/plotly/plotly.rs)), [`charming`](https://crates.io/crates/charming) `0.6.0` ([yuankunzhang/charming](https://github.com/yuankunzhang/charming)) — none ggplot-shaped
  - Existing Rust kind: `rust-native` (the Rust plot libs are independent designs, not ggplot ports)
  - Existing non-C alternatives: `Vega-Lite` (JSON-based, language-agnostic)
  - Parallelism: rayon-able
  - SIMD: none
  - Quadrant: ④
  - GPU-amenable: no — vector plotting
  - Upstream license: `MIT`
  - Priority: `P2`
  - Layer: `adopt` (use Rust plot libs directly; defer grammar-of-graphics)
  - Consumes primitives: —
  - Notes: A grammar-of-graphics Rust plotting library is a worthwhile open project but a big one. Most rsomics workflows can emit Vega-Lite JSON or use `plotly`; defer the grammar-of-graphics dream.

- [ ] **`Plotly` / `plotly.js`** — declarative interactive charts.
  - Reference impl: `JavaScript` · [plotly/plotly.js](https://github.com/plotly/plotly.js) · `MIT`
  - Existing Rust: [`plotly`](https://crates.io/crates/plotly) `0.14.1` — official Rust bindings for plotly.js
  - Existing Rust kind: `pure-port` (faithful Rust bindings emitting plotly.js JSON)
  - Existing non-C alternatives: —
  - Parallelism: browser-side
  - SIMD: none
  - Quadrant: ④
  - GPU-amenable: no — vector charts
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Use `plotly` for interactive HTML charts from Rust pipelines. Mature, well-maintained.

- [ ] **`ETE3` / `ETE4`** — Python tree visualization.
  - Reference impl: `Python` · [etetoolkit/ete](https://github.com/etetoolkit/ete) · `GPL-3`
  - Existing Rust: none verified for visualization; Newick parsers exist in the `phylo` ecosystem
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `iTOL` (web), `FigTree` (Java)
  - Parallelism: Python
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — tree viz
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Tree-viz is GUI work. Provide stable Newick/NEXUS IO and SVG renderers from Rust; leave interactive viewing to ETE / iTOL / FigTree.

- [ ] **`Bandage`** — assembly-graph viewer (Qt desktop).
  - Reference impl: `C++` (Qt) · [rrwick/Bandage](https://github.com/rrwick/Bandage) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `BandageNG` (community fork)
  - Parallelism: Qt threading
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — GUI
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Standalone GUI; not a Rust rewrite target. Make sure our assembler outputs valid GFA so Bandage can read them.

- [ ] **gnomAD viewer / browser**.
  - Reference impl: `React/JavaScript` · [broadinstitute/gnomad-browser](https://github.com/broadinstitute/gnomad-browser) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: browser
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — web app
  - Upstream license: `MIT`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Domain-specific web app; not a port target. Rust pipelines that produce gnomAD-compatible VCF outputs are enough.
