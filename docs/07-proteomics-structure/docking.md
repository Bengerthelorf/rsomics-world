# Molecular docking

> Small-molecule docking, scoring, and binding-pose prediction, plus the
> emerging class of diffusion / structure-prediction-based docking tools.

## Scope

Includes: classical empirical / physics-based docking (AutoDock Vina, Smina,
rDock, Glide), CNN-rescored docking (GNINA), and structure-prediction-style
diffusion docking (DiffDock, RoseTTAFold-AllAtom). Excludes: free-energy
methods (FEP, TI) — out of scope at this stage.

## Design notes

- AutoDock Vina is the open-source workhorse: BSD-3, ~10k LOC C++, well-
  documented scoring function, and the basis for Smina, GNINA, and many
  forks. A pure-Rust Vina is a tractable, high-leverage target.
- GNINA replaces Vina's empirical scoring with a CNN. Inference path is
  small — could be packaged as `candle`/`burn`. Real bottleneck is the
  training data and pretrained weights, not the code.
- DiffDock is a graph-diffusion model for blind docking. Pure-Rust
  inference is feasible via `candle` but the model is non-trivial
  (equivariant message passing). Phase 4+ research project.
- Commercial tools (Glide / Schrödinger Suite) are closed source and
  industry-dominant in pharma; no port path.
- Receptor and ligand prep (PDBQT generation, hydrogens, partial charges)
  is half the practical pain. AutoDockTools, MGLTools, OpenBabel — none
  Rust-native. `openbabel-rs` exists as FFI. Worth a focused crate
  `rsomics-mol-prep` that wraps OpenBabel via FFI.
- License watch: AutoDock Vina **Apache-2.0**, Smina **Apache-2.0 + GPL-2**,
  GNINA **Apache-2.0** (depends on libmolgrid), rDock **LGPL-3**, DiffDock
  **MIT**, RoseTTAFold-AllAtom **MIT** code + non-commercial weights,
  Glide **proprietary commercial**.

## TODO

- [ ] **`AutoDock Vina`** — empirical-scoring small-molecule docking (the open default).
  - Reference impl: `C++` · [ccsb-scripps/AutoDock-Vina](https://github.com/ccsb-scripps/AutoDock-Vina) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `Smina`, `GNINA` (forks)
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — scoring function vectorises; conformer sampling latency-bound
  - Upstream license: `Apache-2.0`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-vina`)
  - Consumes primitives: `pdbtbx`, future `rsomics-mol-prep` (OpenBabel FFI), `nalgebra`, `rayon`
  - Notes: Permissive license, small codebase, published scoring function. Pure-Rust port is a clean, well-scoped project. Inner loops benefit from SIMD.

- [ ] **`Smina`** — fork of Vina with improved scoring-function development support.
  - Reference impl: `C++` · [mwojcikowski/smina](https://github.com/mwojcikowski/smina) · `Apache-2.0 + GPL-2`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `AutoDock Vina`, `GNINA`
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — same as Vina
  - Upstream license: `Apache-2.0 + GPL-2`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-vina` (scoring-function-development mode)
  - Consumes primitives: same as Vina
  - Notes: Once `rsomics-vina` exists, Smina's scoring-function-development affordances slot in as features.

- [ ] **`GNINA`** — Vina + CNN-based scoring.
  - Reference impl: `C++` (libtorch / CUDA) · [gnina/gnina](https://github.com/gnina/gnina) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `DiffDock` (different paradigm)
  - Parallelism: libtorch CUDA
  - SIMD: libtorch kernels
  - Quadrant: —
  - GPU-amenable: yes — CNN rescoring is dense DL
  - Upstream license: `Apache-2.0`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-vina` (CNN-rescoring mode)
  - Consumes primitives: `rsomics-vina`, `candle` or `burn`
  - Notes: Built on Smina + libtorch. The CNN is small. Rust port = `rsomics-vina` core + `candle`/`burn` for the rescoring head. Phase 4+.

- [ ] **`rDock`** — high-throughput virtual-screening docking.
  - Reference impl: `C++` · [rxdock/rxdock](https://gitlab.com/rxdock/rxdock) (community maintenance, GitLab-hosted) · `LGPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `AutoDock Vina`
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — virtual screening trivially parallel
  - Upstream license: `LGPL-3`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Specializes in virtual screening throughput and protein-RNA docking. Niche. GitLab-hosted (gh aliveness N/A).

- [ ] **`DiffDock`** — diffusion-model blind docking.
  - Reference impl: `Python` (PyTorch + PyG) · [gcorso/DiffDock](https://github.com/gcorso/DiffDock) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: PyTorch GPU
  - SIMD: PyTorch kernels
  - Quadrant: —
  - GPU-amenable: yes — equivariant GNN diffusion, dense DL
  - Upstream license: `MIT`
  - Priority: `P2`
  - Layer: `B` (tool — `rsomics-diffdock`)
  - Consumes primitives: `candle` or `burn`, `pdbtbx`, future `rsomics-mol-prep`
  - Notes: Equivariant GNN diffusion. Pure-Rust inference is a candle/burn research project; PyG-style ops not yet first-class in `candle`. Phase 4+.

- [ ] **`RoseTTAFold-AllAtom`** — biomolecular complex prediction (covers docking-like outputs).
  - Reference impl: `Python` (PyTorch) · [baker-laboratory/RoseTTAFold-All-Atom](https://github.com/baker-laboratory/RoseTTAFold-All-Atom) · `MIT` code + non-commercial weights
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `AlphaFold3`, `Boltz-1`, `Chai-1`
  - Parallelism: PyTorch GPU
  - SIMD: PyTorch kernels
  - Quadrant: —
  - GPU-amenable: yes — same family as AF3
  - Upstream license: `MIT` (code) + non-commercial weights
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-boltz` (AF3-class umbrella from [`structure-prediction.md`](structure-prediction.md))
  - Consumes primitives: see Boltz-1 entry
  - Notes: Overlaps heavily with AF3/Boltz/Chai (covered in `structure-prediction.md`). Don't duplicate the port effort. **Cross-reference only — canonical entry is in `structure-prediction.md`.**

- [ ] **`Glide`** (Schrödinger) — commercial industry-standard docking.
  - Reference impl: closed-source · schrodinger.com · commercial
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: closed
  - SIMD: closed
  - Quadrant: —
  - GPU-amenable: unknown — closed
  - Upstream license: commercial
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: No port path; closed source. Listed for completeness.
