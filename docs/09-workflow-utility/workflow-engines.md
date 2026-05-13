# Workflow engines

> Pipeline / workflow orchestration systems for bioinformatics: scheduling,
> dependency resolution, retry, and reproducibility across HPC, cloud, and
> local execution.

## Scope

Includes: rule-based engines (Snakemake), dataflow engines (Nextflow),
WDL-based engines (Cromwell), GUI / web workflows (Galaxy), CWL
implementations (cwltool, Toil), data-lineage engines (Pachyderm, Reflow),
SDK-based managed platforms (Latch), and any nascent Rust workflow
engines. Excludes: low-level job schedulers (Slurm, Kubernetes) and
environment management (see [containers](containers.md)).

## Design notes

- The workflow-engine landscape is mature and entrenched. **Snakemake and
  Nextflow together account for the vast majority of new bioinformatics
  pipelines**.
- Rust's structural advantages (single static binary, type safety,
  fearless concurrency) would be valuable in a workflow engine — but
  none of the existing Rust attempts have community traction.
- A more realistic Rust play is to ensure **first-class integration** with
  Snakemake and Nextflow.
- License watch: Snakemake **MIT**, Nextflow **Apache-2.0**, Cromwell
  **BSD-3**, Galaxy **AFL-3**, cwltool **Apache-2.0**, Toil **Apache-2.0**,
  Pachyderm **Apache-2.0**, Reflow **Apache-2.0**, Latch SDK **Apache-2.0**.

## TODO

- [ ] **`Snakemake`** — Python-based rule + DAG workflow engine.
  - Reference impl: `Python` · [snakemake/snakemake](https://github.com/snakemake/snakemake) · `MIT`
  - Existing Rust: none verified (Snakemake itself supports embedded Rust scripts since v8)
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing + various backends (slurm/k8s)
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — orchestration layer
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: —
  - Consumes primitives: rsomics-* binaries are invoked as rules
  - Notes: Adopt as the recommended user-facing engine. Goal is a `rsomics-snakemake-wrappers` set: ready-made rules for common `rsomics-*` binaries with sensible defaults. Don't build a competitor.

- [ ] **`Nextflow`** — Groovy/JVM dataflow workflow engine.
  - Reference impl: `Groovy/Java` · [nextflow-io/nextflow](https://github.com/nextflow-io/nextflow) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: JVM + various executors
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — orchestration
  - Upstream license: `Apache-2.0`
  - Priority: `P1`
  - Layer: —
  - Consumes primitives: rsomics-* binaries as nf-core modules
  - Notes: Same shape as Snakemake. Provide `rsomics-*` modules for the `nf-core` ecosystem.

- [ ] **`WDL / Cromwell`** — Workflow Description Language + Broad's engine.
  - Reference impl: `Scala` (Cromwell) + WDL spec · [broadinstitute/cromwell](https://github.com/broadinstitute/cromwell) · `BSD-3-Clause`
  - Existing Rust: [`stjude-rust-labs/wdl`](https://github.com/stjude-rust-labs/wdl) — Rust WDL parser, **upstream repo archived 2025-12**
  - Existing Rust kind: `none` (active)
  - Existing non-C alternatives: `miniwdl` (Python)
  - Parallelism: JVM
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — orchestration
  - Upstream license: `BSD-3-Clause`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-workflow` (if a Rust workflow engine ships, WDL parsing mode lives here)
  - Consumes primitives: future `rsomics-wdl-parser`
  - Notes: WDL is mostly used inside Broad / GATK workflows. The stjude-rust-labs WDL parser is now archived; would need fresh effort. Focus on parsing + validation, not execution.

- [ ] **`Galaxy`** — web-based GUI workflow platform.
  - Reference impl: `Python` · [galaxyproject/galaxy](https://github.com/galaxyproject/galaxy) · `AFL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python + various backends
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — web UI
  - Upstream license: `AFL-3`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: rsomics-* binaries as Galaxy ToolShed wrappers
  - Notes: Provide Galaxy tool wrappers for `rsomics-*` binaries. Don't reimplement Galaxy.

- [ ] **`CWL / cwltool`** — Common Workflow Language reference implementation.
  - Reference impl: `Python` · [common-workflow-language/cwltool](https://github.com/common-workflow-language/cwltool) · `Apache-2.0`
  - Existing Rust: [`onnovalkering/cwl`](https://github.com/onnovalkering/cwl) — CWL object model, **upstream repo archived 2025-11**
  - Existing Rust kind: `none` (active)
  - Existing non-C alternatives: `Toil` (Python), `arvados/cwl-runner` (Python)
  - Parallelism: Python multiprocessing
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — orchestration
  - Upstream license: `Apache-2.0`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-workflow` (CWL parsing + execution mode)
  - Consumes primitives: future `rsomics-cwl-parser`, future `rsomics-workflow`
  - Notes: The Rust CWL object-model crate is archived; would need fresh effort. A pure-Rust CWL runner would be useful for static-binary deployment in pharma/clinical settings. Phase 5+.

- [ ] **`Toil`** — Python implementation of CWL/WDL + Python workflows.
  - Reference impl: `Python` · [DataBiosphere/toil](https://github.com/DataBiosphere/toil) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — orchestration
  - Upstream license: `Apache-2.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Niche in pharma + cloud. Same interop-not-rewrite advice.

- [ ] **`Pachyderm`** — data-versioned containerized pipelines.
  - Reference impl: `Go` · [pachyderm/pachyderm](https://github.com/pachyderm/pachyderm) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Go goroutines
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — orchestration
  - Upstream license: `Apache-2.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Niche in bioinformatics; primarily used in industry. Skip.

- [ ] **`Reflow`** — incrementally-cached AWS-native workflows (Grail).
  - Reference impl: `Go` · [grailbio/reflow](https://github.com/grailbio/reflow) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Go goroutines
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — orchestration
  - Upstream license: `Apache-2.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Largely Grail-internal. Skip.

- [ ] **`Latch SDK`** — managed-platform Python SDK for workflows.
  - Reference impl: `Python` · [latchbio/latch](https://github.com/latchbio/latch) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python + managed platform
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — vendor SDK
  - Upstream license: `Apache-2.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Vendor SDK. Provide Latch wrappers for `rsomics-*` if there's user demand.

- [ ] **Rust-native workflow engines** — exploratory.
  - Reference impl: various early-stage projects · `none verified` as production-ready · mixed licenses
  - Existing Rust: no bioinformatics-targeted Rust engine has traction; generic crates exist but are domain-agnostic
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: depends on crate
  - SIMD: depends on crate
  - Quadrant: —
  - GPU-amenable: no — orchestration
  - Upstream license: mixed
  - Priority: `P2`
  - Layer: `B` (future tool — `rsomics-workflow`)
  - Consumes primitives: `tokio`, `petgraph` (DAG), `polars` (provenance), future `rsomics-wdl-parser` / `rsomics-cwl-parser`
  - Notes: No Rust workflow engine has bioinformatics adoption as of 2026-05. A Rust-native lightweight engine is interesting once `rsomics-*` CLIs stabilize, but Phase 5+. Track but don't commit yet.
