# Workflow engines

> Pipeline / workflow orchestration systems for bioinformatics: scheduling,
> dependency resolution, retry, and reproducibility across HPC, cloud, and
> local execution.

## Scope

Includes: rule-based engines (Snakemake), dataflow engines (Nextflow),
WDL-based engines (Cromwell), GUI / web workflows (Galaxy), CWL
implementations (cwltool, Toil), data-lineage engines (Pachyderm, Reflow),
SDK-based managed platforms (Latch), and any nascent Rust workflow
engines (snakegrass, cwl-rs). Excludes: low-level job schedulers (Slurm,
Kubernetes) and environment management (see [containers](containers.md)).

## Design notes

- The workflow-engine landscape is mature and entrenched. **Snakemake and
  Nextflow together account for the vast majority of new bioinformatics
  pipelines**, with WDL/Cromwell dominant inside the Broad / GATK
  ecosystem. Displacing any of these is a multi-year community effort,
  not a code project.
- Rust's structural advantages (single static binary, type safety,
  fearless concurrency) would be valuable in a workflow engine — but
  none of the existing Rust attempts (snakegrass, cwl-rs object model,
  `workflow-engine` crate) have community traction in bioinformatics.
- A more realistic Rust play is to ensure **first-class integration** with
  Snakemake and Nextflow: Snakemake 8+ already supports embedded Rust
  scripts; Nextflow can call any binary, so `rsomics-*` CLIs slot in
  naturally. This means cargo-built static binaries with stable CLI
  contracts are more valuable than yet another engine.
- For the long term, a Rust-native engine has a niche if it solves the
  "WDL/Nextflow are slow to launch hundreds of thousands of tiny tasks
  on local hardware" problem — Snakemake feels this too. A lightweight,
  embedded engine designed for in-process dispatch from a Rust pipeline
  could be valuable. Phase 5+.
- License watch: Snakemake **MIT**, Nextflow **Apache-2.0**, Cromwell
  **BSD-3**, Galaxy **AFL-3**, cwltool **Apache-2.0**, Toil **Apache-2.0**,
  Pachyderm **Apache-2.0**, Reflow **Apache-2.0**, Latch SDK **Apache-2.0**.

## TODO

- [ ] **`Snakemake`** — Python-based rule + DAG workflow engine.
  - Reference impl: `Python` · [snakemake/snakemake](https://github.com/snakemake/snakemake) · `MIT`
  - Existing Rust: none verified (Snakemake itself supports embedded Rust scripts since v8)
  - Existing non-C alternatives: —
  - Priority: `P1` (interop)
  - Notes: Adopt as the recommended user-facing engine. Goal is a
    `rsomics-snakemake-wrappers` set: ready-made rules for common
    `rsomics-*` binaries with sensible defaults. Don't build a competitor.

- [ ] **`Nextflow`** — Groovy/JVM dataflow workflow engine.
  - Reference impl: `Groovy/Java` · [nextflow-io/nextflow](https://github.com/nextflow-io/nextflow) · `Apache-2.0`
  - Existing Rust: none
  - Existing non-C alternatives: —
  - Priority: `P1` (interop)
  - Notes: Same shape as Snakemake. Provide `rsomics-*` modules for the
    `nf-core` ecosystem; aim for upstream-acceptable PRs that add Rust-
    based tools as nf-core modules.

- [ ] **`WDL / Cromwell`** — Workflow Description Language + Broad's engine.
  - Reference impl: `Scala` (Cromwell) + WDL spec · [broadinstitute/cromwell](https://github.com/broadinstitute/cromwell) · `BSD-3-Clause`
  - Existing Rust: [`stjude-rust-labs/sprocket` / `wdl`](https://github.com/stjude-rust-labs/wdl) — Rust WDL parser + linter (verify; mentioned in literature). If unverified, mark as `none verified`.
  - Existing non-C alternatives: `miniwdl` (Python)
  - Priority: `P2`
  - Notes: WDL is mostly used inside Broad / GATK workflows. A Rust WDL
    *parser* is realistic (compare to the cwl-rs object model); a Rust WDL
    *executor* duplicates Cromwell. Focus on parsing + validation, not
    execution.

- [ ] **`Galaxy`** — web-based GUI workflow platform.
  - Reference impl: `Python` · [galaxyproject/galaxy](https://github.com/galaxyproject/galaxy) · `AFL-3`
  - Existing Rust: none
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Provide Galaxy tool wrappers for `rsomics-*` binaries so they
    appear in the Galaxy ToolShed. Don't reimplement Galaxy.

- [ ] **`CWL / cwltool`** — Common Workflow Language reference implementation.
  - Reference impl: `Python` · [common-workflow-language/cwltool](https://github.com/common-workflow-language/cwltool) · `Apache-2.0`
  - Existing Rust: [`onnovalkering/cwl`](https://github.com/onnovalkering/cwl) — CWL object model in Rust (parser, no executor)
  - Existing non-C alternatives: `Toil` (Python), `arvados/cwl-runner` (Python)
  - Priority: `P2`
  - Notes: The Rust CWL object-model crate is a foundation, not an engine.
    A pure-Rust CWL runner would be useful for static-binary deployment
    in pharma/clinical settings. Phase 5+.

- [ ] **`Toil`** — Python implementation of CWL/WDL + Python workflows.
  - Reference impl: `Python` · [DataBiosphere/toil](https://github.com/DataBiosphere/toil) · `Apache-2.0`
  - Existing Rust: none
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Niche in pharma + cloud. Same interop-not-rewrite advice.

- [ ] **`Pachyderm`** — data-versioned containerized pipelines.
  - Reference impl: `Go` · [pachyderm/pachyderm](https://github.com/pachyderm/pachyderm) · `Apache-2.0`
  - Existing Rust: none
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Niche in bioinformatics; primarily used in industry. Skip.

- [ ] **`Reflow`** — incrementally-cached AWS-native workflows (Grail).
  - Reference impl: `Go` · [grailbio/reflow](https://github.com/grailbio/reflow) · `Apache-2.0`
  - Existing Rust: none
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Largely Grail-internal. Skip.

- [ ] **`Latch SDK`** — managed-platform Python SDK for workflows.
  - Reference impl: `Python` · [latchbio/latch](https://github.com/latchbio/latch) · `Apache-2.0`
  - Existing Rust: none
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Vendor SDK. Provide Latch wrappers for `rsomics-*` if there's
    user demand; otherwise skip.

- [ ] **Rust-native workflow engines** — exploratory.
  - Reference impl: various early-stage projects · `none verified` as production-ready · mixed licenses
  - Existing Rust: [`workflow-engine` on docs.rs](https://docs.rs/workflow-engine/latest/workflow_engine/) — generic engine, not bioinformatics-specific; `onnovalkering/cwl` — CWL object model only
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: No Rust workflow engine has bioinformatics adoption as of
    2026-05. A Rust-native lightweight engine is interesting once
    `rsomics-*` CLIs stabilize, but Phase 5+. Track but don't commit yet.
