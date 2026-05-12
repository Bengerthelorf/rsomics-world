# 09 — Workflow & utility

> Cross-cutting concerns: workflow / pipeline orchestration, reproducible
> environments and containers, and data visualization.

## Sub-areas

| File | Scope |
|------|-------|
| [`workflow-engines.md`](workflow-engines.md) | Snakemake, Nextflow, WDL/Cromwell, Galaxy, CWL/cwltool, Toil, Pachyderm, Reflow, Latch SDK, and Rust experiments (snakegrass, `cwl-rs`). |
| [`containers.md`](containers.md) | Package and environment management: Conda, Mamba, Pixi (Rust!), Bioconda, Singularity/Apptainer, Docker, Wave, BioContainers, environment-modules. |
| [`data-viz.md`](data-viz.md) | Genome browsers, plotting libraries, and assembly-graph viewers: IGV, JBrowse2, UCSC Browser, Gosling, deepTools, ggplot2 (R), Plotly, ETE3, Bandage, gnomAD viewer. |

## Design notes

- This module is **horizontal** — every other module's tools eventually have
  to be wrapped in a workflow, packaged in a container, and visualized.
  Getting workflow + env management right multiplies the leverage of every
  rewrite in modules 01-08.
- Pixi is the standout Rust success story in this space: a fast, reproducible,
  conda-compatible package manager written entirely in Rust. Adopt as the
  recommended environment manager for `rsomics-world` itself.
- Workflow engines remain Python (Snakemake), Groovy/JVM (Nextflow), or
  Java (Cromwell). Rust-native workflow engines exist but none are mature
  or widely adopted in bioinformatics as of 2026. This is a real
  ecosystem gap with no obvious winner yet.
- Visualization is mostly JavaScript (IGV.js, JBrowse2, Gosling) or
  desktop GUIs (PyMOL, ChimeraX, Bandage) — Rust has limited play except
  for headless rendering backends and bulk plot generation. WebAssembly
  is the realistic Rust-in-the-browser path.
- The realistic Rust deliverables in this module are: (a) excellent
  `cargo`-based packaging that complements pixi; (b) Rust libraries for
  reading/writing common viz file formats (BigWig, BigBed — there is
  `bigtools` already in Rust); (c) a polished CLI experience that
  reduces the need for workflow engines for simple cases.
- License watch: Snakemake **MIT**, Nextflow **Apache-2.0**, Cromwell
  **BSD-3**, Galaxy **AFL-3**, cwltool **Apache-2.0**, Toil **Apache-2.0**,
  Pixi **BSD-3**, Conda **BSD-3**, Mamba **BSD-3**, Apptainer **BSD-3**,
  Docker **Apache-2.0**, IGV **MIT**, JBrowse2 **Apache-2.0**.
