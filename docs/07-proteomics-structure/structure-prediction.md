# Protein structure prediction

> AI-based 3D structure prediction of proteins and biomolecular complexes
> from sequence (and, increasingly, ligand/multimer context).

## Scope

Includes: MSA-based predictors (AlphaFold2, OpenFold, RoseTTAFold,
ColabFold), MSA-free / language-model predictors (ESMFold, OmegaFold), and
the next-generation diffusion-style multi-molecule predictors (AlphaFold3,
Boltz-1, Chai-1, RoseTTAFold-AllAtom). Excludes: structure visualization
and analysis (see [structure-analysis](structure-analysis.md)) and
docking (see [docking](docking.md)).

## Design notes

- This entire sub-area is PyTorch-bound. Pure-Rust implementations of any
  of these models are research-grade only and bottlenecked on (a) reading
  the publicly released weights, (b) reproducing the inference graph, and
  (c) handling the MSA-generation upstream. The **bottleneck is the
  weights and the science, not the code language**.
- `candle` (HuggingFace) and `burn` are both production-quality enough to
  load and run these architectures for inference. `candle-transformers`
  already ships ESM-2. AlphaFold-class inference is feasible but nobody
  has shipped a polished pure-Rust binary yet.
- MSA generation (HHblits, JackHMMER) is the upstream cost. ColabFold uses
  MMseqs2 against a remote service to bypass this. A Rust MSA-generation
  toolchain (Rust `MMseqs2`-equivalent) would help every tool in this
  list and is independently valuable; see also classification.md.
- For Phase 4+, the realistic deliverables are: (a) a `candle`-based
  ESMFold inference binary; (b) a Rust ColabFold-style MSA orchestrator;
  (c) wrappers that present a unified CLI over AF2/AF3/Boltz when the
  user has GPUs. Training is out of scope.
- License watch: AlphaFold2 code **Apache-2.0** (weights CC BY 4.0),
  AlphaFold3 **CC BY-NC-SA-4.0** (weights separate, non-commercial),
  ColabFold **MIT**, ESMFold **MIT** (Meta), RoseTTAFold **MIT** code +
  Rosetta-DL non-commercial weights, OmegaFold check repo, OpenFold
  **Apache-2.0**, OpenFold-3 Apache-2.0, Boltz-1 **MIT**, Chai-1 **Apache-2.0**.

## TODO

- [ ] **`AlphaFold2`** — DeepMind's MSA-based protein structure predictor.
  - Reference impl: `Python` (JAX) · [google-deepmind/alphafold](https://github.com/google-deepmind/alphafold) · `Apache-2.0` (code) / `CC BY 4.0` (weights)
  - Existing Rust: none verified for full inference
  - Existing non-C alternatives: `OpenFold` (PyTorch reimplementation), `ColabFold` (deployment)
  - Priority: `P1`
  - Notes: Pure-Rust inference is a multi-month research effort with
    `candle`/`burn`. The community is more likely to use AF2 via
    ColabFold (which is `Python` + `MMseqs2`); Rust value is in providing
    a fast MSA-generation backend, not the network.

- [ ] **`AlphaFold3`** — DeepMind's diffusion-based multi-molecule predictor.
  - Reference impl: `Python` (JAX) · [google-deepmind/alphafold3](https://github.com/google-deepmind/alphafold3) · `CC BY-NC-SA-4.0` (code) + AF3 model-parameters ToU (non-commercial, gated)
  - Existing Rust: none verified
  - Existing non-C alternatives: `OpenFold-3` (Apache-2.0 reproduction), `Boltz-1`, `Chai-1`
  - Priority: `P2` (weights are gated; commercial pure-Rust users will use Boltz/Chai instead)
  - Notes: Non-commercial license + gated weights make AF3 itself a poor
    Rust port target. For commercial users, target `Boltz-1` or `Chai-1`
    (Apache/MIT) instead.

- [ ] **`ColabFold`** — community deployment of AF2 with MMseqs2 MSAs.
  - Reference impl: `Python` · [sokrypton/ColabFold](https://github.com/sokrypton/ColabFold) · `MIT` (code) + Apache-2.0 (AF2 inherited)
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: ColabFold's contribution is the MMseqs2 remote-MSA service +
    Colab UX, not the inference. A Rust analogue is "a `mmseqs2-rs` MSA
    pipeline that emits AF2-compatible feature dicts". Useful glue.

- [ ] **`ESMFold`** — Meta's MSA-free language-model-based predictor.
  - Reference impl: `Python` (PyTorch) · [facebookresearch/esm](https://github.com/facebookresearch/esm) · `MIT`
  - Existing Rust: partial — `candle-transformers` ships ESM-2 weights inference; ESMFold-specific head not yet packaged
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Best candidate for the first pure-Rust structure-prediction
    binary. No MSA dependency, single-sequence forward pass. Weights and
    architecture are open. Extend `candle-transformers` to add the
    ESMFold folding head; this is a tractable Phase-4 project.

- [ ] **`RoseTTAFold`** / `RoseTTAFold2` — Baker-lab's three-track architecture.
  - Reference impl: `Python` (PyTorch) · [RosettaCommons/RoseTTAFold](https://github.com/RosettaCommons/RoseTTAFold) · `MIT` (code) + Rosetta-DL non-commercial weights
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Strong scientifically but weights are non-commercial; users
    with commercial needs avoid it. AF2/ESMFold/Boltz cover the rest.

- [ ] **`OmegaFold`** — MSA-free language-model predictor (HelixonAI).
  - Reference impl: `Python` (PyTorch) · [HeliXonProtein/OmegaFold](https://github.com/HeliXonProtein/OmegaFold) · check repo (Apache-2.0 reported in literature)
  - Existing Rust: none verified
  - Existing non-C alternatives: `ESMFold` (Meta, similar idea)
  - Priority: `P2`
  - Notes: Same shape as ESMFold but less ecosystem traction. Skip unless
    a specific benchmark needs it.

- [ ] **`OpenFold`** — community PyTorch reproduction of AlphaFold2.
  - Reference impl: `Python` (PyTorch) · [aqlaboratory/openfold](https://github.com/aqlaboratory/openfold) · `Apache-2.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Easier port target than AlphaFold itself — clean PyTorch with
    Apache-2.0 weights. If someone is going to build pure-Rust AF2
    inference, start from OpenFold's architecture, not DeepMind's JAX.

- [ ] **`Boltz-1`** — open AF3-class diffusion predictor.
  - Reference impl: `Python` (PyTorch) · [jwohlwend/boltz](https://github.com/jwohlwend/boltz) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: `Chai-1` (similar role)
  - Priority: `P1`
  - Notes: The de-facto open AF3-class model in 2025-2026 — MIT, weights
    open, no commercial restriction. Long-term most important target for a
    commercial-friendly pure-Rust AF3-class inference binary.

- [ ] **`Chai-1`** — Chai Discovery's open AF3-class predictor.
  - Reference impl: `Python` (PyTorch) · [chaidiscovery/chai-lab](https://github.com/chaidiscovery/chai-lab) · `Apache-2.0` (code) (weights subject to terms)
  - Existing Rust: none verified
  - Existing non-C alternatives: `Boltz-1`
  - Priority: `P2`
  - Notes: Choose either Boltz-1 or Chai-1 as the open-AF3-class target;
    we don't need both.
