# Containers and environment management

> Reproducible software environments: package managers, container runtimes,
> bioinformatics package channels, and HPC module systems.

## Scope

Includes: conda-family package managers (Conda, Mamba, Pixi), bioinformatics
package channels (Bioconda, BioContainers), container runtimes (Docker,
Singularity/Apptainer), container-on-demand builders (Wave from Nextflow),
and HPC environment-modules. Excludes: cluster orchestration / Kubernetes.

`rsomics-world` itself recommends **pixi + a `pixi.toml` per analysis**
as the default user-facing environment manager, with `cargo` for any
Rust code. We co-distribute binaries via Bioconda so they coexist with
the existing ecosystem.

## Design notes

- Pixi is the rare case of a Rust rewrite that won. Conda-compatible,
  much faster than `conda`/`mamba` for solving, has cargo-style lockfiles,
  and is genuinely cross-platform. Adopt as our recommended user-facing
  env manager.
- Conda and Mamba are not replacement targets — they are baseline.
- Apptainer is the de-facto HPC container runtime. No Rust play.
- Docker is the container runtime. No Rust play (`youki` is a Rust OCI
  runtime, but not what bioinformatics users interact with).
- License watch: Pixi **BSD-3**, Conda **BSD-3**, Mamba **BSD-3**,
  Apptainer **BSD-3**, Docker (Engine) **Apache-2.0**, Wave **Apache-2.0**.

## TODO

- [x] **`Pixi`** — fast Rust-native conda-compatible package manager.
  - Reference impl: `Rust` · [prefix-dev/pixi](https://github.com/prefix-dev/pixi) · `BSD-3-Clause`
  - Existing Rust: [`pixi`](https://crates.io/crates/pixi) `0.15.2`
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: —
  - Parallelism: tokio
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: no — package management
  - Upstream license: `BSD-3-Clause`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Recommended env manager for users of `rsomics-*`. Ship a `pixi.toml` template in this repo and reference it from each crate's README.

- [ ] **`Conda`** — original Python-centric package manager.
  - Reference impl: `Python` · [conda/conda](https://github.com/conda/conda) · `BSD-3-Clause`
  - Existing Rust: none direct (Pixi reads conda packages but is not a port)
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `Mamba`, `Pixi`
  - Parallelism: Python multiprocessing
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — package management
  - Upstream license: `BSD-3-Clause`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Baseline; not a rewrite target. Ensure compatibility.

- [ ] **`Mamba` / `micromamba`** — C++ reimplementation of conda's solver.
  - Reference impl: `C++` (libsolv) · [mamba-org/mamba](https://github.com/mamba-org/mamba) · `BSD-3-Clause`
  - Existing Rust: [`rattler`](https://crates.io/crates/rattler) `0.43.0` ([conda/rattler](https://github.com/conda/rattler)) — Rust conda toolkit underpinning Pixi
  - Existing Rust kind: `rust-native` (rattler reimplements the conda package handling in pure Rust, not a code-port of mamba)
  - Existing non-C alternatives: `Pixi`
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: ①
  - GPU-amenable: no — package management
  - Upstream license: `BSD-3-Clause`
  - Priority: `P2`
  - Layer: `adopt` (rattler as library; pixi as CLI)
  - Consumes primitives: —
  - Notes: Rust users should use `rattler` directly (as a library) or `pixi` (as a CLI). `mamba` itself is fine as-is; not a rewrite target.

- [x] **`rattler`** — Rust conda-package toolkit (underlies Pixi).
  - Reference impl: `Rust` · [conda/rattler](https://github.com/conda/rattler) · `BSD-3-Clause`
  - Existing Rust: [`rattler`](https://crates.io/crates/rattler) `0.43.0`
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: —
  - Parallelism: tokio
  - SIMD: auto-vectorize
  - Quadrant: ①
  - GPU-amenable: no — package management
  - Upstream license: `BSD-3-Clause`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Use this when we need to manipulate conda packages from Rust code (e.g. QA tooling on bioconda recipes).

- [ ] **`Bioconda`** — community bioinformatics conda channel.
  - Reference impl: recipes repo · [bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes) · `MIT`
  - Existing Rust: none (infrastructure)
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: CI-driven
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — recipes repository
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: —
  - Consumes primitives: —
  - Notes: Ship every `rsomics-*` crate as a bioconda recipe at first release. Bioconda's CI handles container builds (BioContainers) for free.

- [ ] **`Singularity` / `Apptainer`** — HPC container runtime.
  - Reference impl: `Go` · [apptainer/apptainer](https://github.com/apptainer/apptainer) · `BSD-3-Clause`
  - Existing Rust: none (cf. `youki` for general OCI runtimes)
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Go goroutines
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — container runtime
  - Upstream license: `BSD-3-Clause`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: We test our containers under Apptainer; no rewrite.

- [ ] **`Docker`** — general container runtime.
  - Reference impl: `Go` · [moby/moby](https://github.com/moby/moby) · `Apache-2.0`
  - Existing Rust: [`youki`](https://github.com/containers/youki) — Rust OCI runtime (general, not bioinformatics)
  - Existing Rust kind: `none` (for Docker engine specifically; youki is a parallel OCI implementation)
  - Existing non-C alternatives: `podman`, `containerd`
  - Parallelism: Go goroutines
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — container runtime
  - Upstream license: `Apache-2.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Not a bioinformatics rewrite target.

- [ ] **`Wave`** — Nextflow's on-demand container builder.
  - Reference impl: `Java` · [seqeralabs/wave](https://github.com/seqeralabs/wave) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: JVM
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — service
  - Upstream license: `Apache-2.0`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Service; we make sure our bioconda recipes are wave-buildable.

- [ ] **`BioContainers`** — community-driven biocontainer registry + recipes.
  - Reference impl: recipes / Docker · [BioContainers](https://github.com/BioContainers) · `Apache-2.0`
  - Existing Rust: none (infrastructure)
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: CI-driven
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — infrastructure
  - Upstream license: `Apache-2.0`
  - Priority: `P1`
  - Layer: —
  - Consumes primitives: —
  - Notes: Same as Bioconda: contribute, don't rebuild.

- [ ] **`environment-modules` / `Lmod`** — classical HPC module system.
  - Reference impl: `Tcl` (classic) / `Lua` (Lmod) · [TACC/Lmod](https://github.com/TACC/Lmod) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: shell
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — module system
  - Upstream license: `MIT`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Legacy HPC infrastructure. Document compatibility; no Rust work.
