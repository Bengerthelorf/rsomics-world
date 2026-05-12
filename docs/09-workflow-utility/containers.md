# Containers and environment management

> Reproducible software environments: package managers, container runtimes,
> bioinformatics package channels, and HPC module systems.

## Scope

Includes: conda-family package managers (Conda, Mamba, Pixi), bioinformatics
package channels (Bioconda, BioContainers), container runtimes (Docker,
Singularity/Apptainer), container-on-demand builders (Wave from Nextflow),
and HPC environment-modules. Excludes: cluster orchestration / Kubernetes,
covered by ops literature.

## Scope notes

`rsomics-world` itself recommends **pixi + a `pixi.toml` per analysis**
as the default user-facing environment manager, with `cargo` for any
Rust code. We co-distribute binaries via Bioconda so they coexist with
the existing ecosystem.

## Design notes

- Pixi is the rare case of a Rust rewrite that won. It is conda-compatible
  (uses the conda ecosystem of packages), much faster than `conda`/`mamba`
  for solving, has cargo-style lockfiles, and is genuinely cross-platform.
  Adopt as our recommended user-facing env manager.
- Conda and Mamba are not replacement targets — they are baseline. We
  don't rewrite them; we ship into the bioconda channel so we coexist
  with them.
- Apptainer (the open-governance fork of Singularity) is the de-facto
  HPC container runtime. No Rust play; we ensure our binaries run
  cleanly in Apptainer images.
- Docker is the container runtime. No Rust play. (`youki` is a Rust OCI
  runtime, but that's not what bioinformatics users interact with.)
- BioContainers and Bioconda are infrastructure projects — we contribute
  recipes, we don't rewrite them. Each `rsomics-*` crate gets a
  conda-forge / bioconda recipe at first release.
- Wave (Seqera / Nextflow) builds containers on demand from conda specs.
  Useful but not something to reimplement; we make sure our packages
  are wave-compatible.
- License watch: Pixi **BSD-3**, Conda **BSD-3**, Mamba **BSD-3**,
  Apptainer **BSD-3**, Docker (Engine) **Apache-2.0**, Wave **Apache-2.0**,
  conda-forge / Bioconda **infrastructure** (recipes individually licensed).

## TODO

- [x] **`Pixi`** — fast Rust-native conda-compatible package manager.
  - Reference impl: `Rust` · [prefix-dev/pixi](https://github.com/prefix-dev/pixi) · `BSD-3-Clause`
  - Existing Rust: same as reference
  - Existing non-C alternatives: —
  - Priority: `P0` (adopt)
  - Notes: Recommended env manager for users of `rsomics-*`. Ship a
    `pixi.toml` template in this repo and reference it from each crate's
    README.

- [ ] **`Conda`** — original Python-centric package manager.
  - Reference impl: `Python` · [conda/conda](https://github.com/conda/conda) · `BSD-3-Clause`
  - Existing Rust: none (Pixi reads conda packages but is not a port)
  - Existing non-C alternatives: `Mamba`, `Pixi`
  - Priority: `P2`
  - Notes: Baseline; not a rewrite target. Ensure compatibility.

- [ ] **`Mamba` / `micromamba`** — C++ reimplementation of conda's solver.
  - Reference impl: `C++` (libsolv) · [mamba-org/mamba](https://github.com/mamba-org/mamba) · `BSD-3-Clause`
  - Existing Rust: `rattler` ([conda/rattler](https://github.com/conda/rattler)) — Rust conda toolkit underpinning Pixi
  - Existing non-C alternatives: `Pixi`
  - Priority: `P2` (Rattler already covers this)
  - Notes: Rust users should use `rattler` directly (as a library) or
    `pixi` (as a CLI). `mamba` itself is fine as-is; not a rewrite target.

- [x] **`rattler`** — Rust conda-package toolkit (underlies Pixi).
  - Reference impl: `Rust` · [conda/rattler](https://github.com/conda/rattler) · `BSD-3-Clause`
  - Existing Rust: same as reference
  - Existing non-C alternatives: —
  - Priority: `P0` (adopt for any programmatic conda handling)
  - Notes: Use this when we need to manipulate conda packages from Rust
    code (e.g., for QA tooling on bioconda recipes).

- [ ] **`Bioconda`** — community bioinformatics conda channel.
  - Reference impl: recipes repo · [bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes) · `MIT` (infrastructure)
  - Existing Rust: none (infrastructure)
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Ship every `rsomics-*` crate as a bioconda recipe at first
    release. Bioconda's CI handles container builds (BioContainers) for free.

- [ ] **`Singularity` / `Apptainer`** — HPC container runtime.
  - Reference impl: `Go` · [apptainer/apptainer](https://github.com/apptainer/apptainer) · `BSD-3-Clause`
  - Existing Rust: none (cf. `youki` for general OCI runtimes)
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: We test our containers under Apptainer; no rewrite.

- [ ] **`Docker`** — general container runtime.
  - Reference impl: `Go` · [moby/moby](https://github.com/moby/moby) · `Apache-2.0`
  - Existing Rust: `youki` ([containers/youki](https://github.com/containers/youki)) — Rust OCI runtime
  - Existing non-C alternatives: `podman`, `containerd`
  - Priority: `P2`
  - Notes: Not a bioinformatics rewrite target.

- [ ] **`Wave`** — Nextflow's on-demand container builder.
  - Reference impl: `Java` · [seqeralabs/wave](https://github.com/seqeralabs/wave) · `Apache-2.0`
  - Existing Rust: none
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Service; we make sure our bioconda recipes are wave-buildable.

- [ ] **`BioContainers`** — community-driven biocontainer registry + recipes.
  - Reference impl: recipes / Docker · [BioContainers](https://github.com/BioContainers) · `Apache-2.0` (infra)
  - Existing Rust: none (infrastructure)
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Same as Bioconda: contribute, don't rebuild.

- [ ] **`environment-modules` / `Lmod`** — classical HPC module system.
  - Reference impl: `Tcl` (classic) / `Lua` (Lmod) · [TACC/Lmod](https://github.com/TACC/Lmod) · `MIT`
  - Existing Rust: none
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Legacy HPC infrastructure. Document compatibility; no
    Rust work.
