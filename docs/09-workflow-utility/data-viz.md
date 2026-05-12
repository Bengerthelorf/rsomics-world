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
  Gosling, igv-reports) and **R** (ggplot2). Rust's play is mostly
  back-end:
  - **Fast bulk renderers** for headless plot generation in pipelines
    (replace deepTools' Python+matplotlib bulk-coverage plotting with
    `plotters` / `plotly-rs` + parallel rendering).
  - **WebAssembly viz components** compiled from Rust for embedding in
    JBrowse2 / IGV.js (still experimental).
  - **Backend services** that serve BigWig / BigBed / CRAM ranges
    efficiently. `bigtools` is already Rust for BigWig/BigBed.
- For tree viz (ETE3, iTOL), Rust has no contender; ETE3 is the Python
  standard. Provide good Newick / NEXUS IO via a canonical Rust tree
  crate so Python tooling can read our outputs.
- The right viz-side rewrites are small, focused crates that the larger
  JS / R ecosystem can call as binaries (e.g., a fast `rsomics-coverage`
  that emits PNG/SVG plots from BAM).
- Bandage is Qt/C++ — a desktop GUI. No Rust play unless we want a Rust
  Tauri-based reimplementation, which is overkill.
- License watch: IGV **MIT**, JBrowse2 **Apache-2.0**, UCSC Genome
  Browser **non-commercial source-available** (UCSC license, free for
  academic), Gosling **MIT**, deepTools **GPL-3**, ggplot2 **MIT**,
  Plotly.js **MIT**, ETE3 **GPL-3**, Bandage **GPL-3**.

## TODO

- [ ] **`IGV` / `IGV.js`** — Integrative Genomics Viewer (desktop + web).
  - Reference impl: `Java` (desktop) + `JavaScript` (web) · [igvteam/igv](https://github.com/igvteam/igv) and [igvteam/igv.js](https://github.com/igvteam/igv.js) · `MIT`
  - Existing Rust: none for viewer itself; backend file-serving via `bigtools`, `noodles`
  - Existing non-C alternatives: `JBrowse2`
  - Priority: `P2`
  - Notes: Don't rewrite the viewer. Provide solid Rust IO so IGV can
    load our outputs (CRAM, BigWig, VCF) over HTTP range requests.

- [ ] **`JBrowse2`** — modern JavaScript / TypeScript genome browser.
  - Reference impl: `TypeScript` · [GMOD/jbrowse-components](https://github.com/GMOD/jbrowse-components) · `Apache-2.0`
  - Existing Rust: none for viewer; WASM-based plugin path exists in JBrowse2
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: WASM-compiled Rust plugins for JBrowse2 are a realistic
    Phase 5+ direction. Until then, focus on backend.

- [ ] **`UCSC Genome Browser`** — classic web genome browser + Kent tools.
  - Reference impl: `C` · [ucscGenomeBrowser/kent](https://github.com/ucscGenomeBrowser/kent) · UCSC source-available (non-commercial; commercial license available)
  - Existing Rust: [`bigtools`](https://github.com/jackh726/bigtools) ★ — pure-Rust BigWig/BigBed reader/writer (replaces parts of `kent`)
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: We do not rebuild the browser. The Rust win is `bigtools` and
    related format readers/writers, which let pipelines avoid the UCSC
    license entirely. Adopt `bigtools`.

- [x] **`bigtools`** — Rust BigWig/BigBed library + CLI.
  - Reference impl: `Rust` · [jackh726/bigtools](https://github.com/jackh726/bigtools) · `MIT`
  - Existing Rust: same as reference
  - Existing non-C alternatives: —
  - Priority: `P0` (adopt)
  - Notes: Foundational. Use it wherever pipelines need BigWig/BigBed.

- [ ] **`Gosling`** — grammar-based interactive genomic viz.
  - Reference impl: `TypeScript` · [gosling-lang/gosling.js](https://github.com/gosling-lang/gosling.js) · `MIT`
  - Existing Rust: none
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: JS only; Rust play is producing Gosling-spec-compatible JSON
    output from analysis pipelines.

- [ ] **`deepTools`** — Python coverage / heatmap plotting from BAM/BigWig.
  - Reference impl: `Python` · [deeptools/deepTools](https://github.com/deeptools/deepTools) · `GPL-3`
  - Existing Rust: none verified (but the pieces exist — `noodles-bam`, `bigtools`, `plotters`)
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Excellent Rust rewrite target. Inner loops (coverage
    aggregation over millions of intervals × many BAMs) are CPU-bound and
    parallelize trivially with `rayon`. Plotting via `plotters` /
    `plotly-rs`. A `rsomics-deeptools`-equivalent would be measurably
    faster and avoid Python deps for pipelines.

- [ ] **`ggplot2`** — R grammar-of-graphics plotting (de facto standard).
  - Reference impl: `R` · [tidyverse/ggplot2](https://github.com/tidyverse/ggplot2) · `MIT`
  - Existing Rust: partial — [`plotters`](https://github.com/plotters-rs/plotters), [`plotly-rs`](https://github.com/plotly/plotly.rs), [`charming`](https://github.com/charming-rs/charming) — none ggplot-shaped
  - Existing non-C alternatives: `Vega-Lite` (JSON-based, language-agnostic)
  - Priority: `P2`
  - Notes: A grammar-of-graphics Rust plotting library is a worthwhile
    open project but a big one (`gglot-rs` attempts exist on GitHub but
    not mature). Most rsomics workflows can emit Vega-Lite JSON or use
    `plotly-rs`; defer the grammar-of-graphics dream.

- [ ] **`Plotly` / `plotly.js`** — declarative interactive charts.
  - Reference impl: `JavaScript` · [plotly/plotly.js](https://github.com/plotly/plotly.js) · `MIT`
  - Existing Rust: [`plotly-rs`](https://github.com/plotly/plotly.rs) ★ — official Rust bindings for plotly.js
  - Existing non-C alternatives: —
  - Priority: `P1` (adopt)
  - Notes: Use `plotly-rs` for interactive HTML charts from Rust pipelines.
    Mature, well-maintained.

- [ ] **`ETE3` / `ETE4`** — Python tree visualization.
  - Reference impl: `Python` · [etetoolkit/ete](https://github.com/etetoolkit/ete) · `GPL-3`
  - Existing Rust: none verified for visualization; Newick parsers exist (`rust-newick`, `phylo` ecosystem)
  - Existing non-C alternatives: `iTOL` (web), `FigTree` (Java)
  - Priority: `P2`
  - Notes: Tree-viz is GUI work. Provide stable Newick/NEXUS IO and SVG
    renderers from Rust; leave interactive viewing to ETE / iTOL / FigTree.

- [ ] **`Bandage`** — assembly-graph viewer (Qt desktop).
  - Reference impl: `C++` (Qt) · [rrwick/Bandage](https://github.com/rrwick/Bandage) · `GPL-3`
  - Existing Rust: none verified
  - Existing non-C alternatives: `BandageNG` (community fork)
  - Priority: `P2`
  - Notes: Standalone GUI; not a Rust rewrite target. Make sure our
    assembler outputs valid GFA so Bandage can read them.

- [ ] **gnomAD viewer / browser**.
  - Reference impl: `React/JavaScript` · [broadinstitute/gnomad-browser](https://github.com/broadinstitute/gnomad-browser) · `MIT`
  - Existing Rust: none
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Domain-specific web app; not a port target. Rust pipelines
    that produce gnomAD-compatible VCF outputs are enough.
