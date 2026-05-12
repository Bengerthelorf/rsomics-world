# 07 — Proteomics & structural biology

> Mass-spectrometry-based peptide/protein identification and quantification,
> protein structure prediction, structure analysis, and molecular docking.

## Sub-areas

| File | Scope |
|------|-------|
| [`mass-spectrometry.md`](mass-spectrometry.md) | DDA + DIA peptide identification, search engines, integrated suites. MaxQuant, MSFragger / FragPipe, Comet, X!Tandem, MSGF+, OpenMS, DIA-NN, Skyline, Spectronaut, Quandenser. |
| [`structure-prediction.md`](structure-prediction.md) | AI-based 3D structure prediction: AlphaFold2/3, ColabFold, ESMFold, RoseTTAFold, OmegaFold, OpenFold, Boltz-1, Chai-1. |
| [`structure-analysis.md`](structure-analysis.md) | Structure file IO, visualization, comparison, secondary-structure and accessibility: PyMOL, ChimeraX, MDAnalysis, BioPython.PDB, ProDy, freesasa, DSSP, Foldseek, TM-align, FATCAT. |
| [`docking.md`](docking.md) | Small-molecule and biomolecular docking: AutoDock Vina, GNINA, Smina, rDock, DiffDock, RoseTTAFold-AllAtom, Glide. |

## Design notes

- This is by far the **least Rust-mature** of all the domains. Mass spec
  tooling is dominated by closed-source Java/C# (MaxQuant, Skyline) and
  vendor binaries; structure prediction is Python+PyTorch with multi-GB
  weights; docking is C++. We will be doing more wrapping and less
  rewriting here than in any other module.
- The ML tools (AlphaFold family, ESMFold, RoseTTAFold, DiffDock) all live
  on PyTorch. Pure-Rust inference is feasible via `candle` or `burn` — the
  Rust core team and HuggingFace already ship `candle-transformers` with
  many of these architectures. **Bottleneck is models + weights, not code.**
  Phase-staged plan: Phase 4+ for any pure-Rust inference; Phase 1
  contributions are limited to data-prep, MSA generation, and post-
  processing utilities.
- mzML/mzXML/MGF parsing is the equivalent of "noodles" for proteomics.
  `mzdata` (Rust) exists and is the obvious foundation. Build on top of it.
- Foldseek is already a hand-tuned C++/SIMD tool by Steineggerlab; it's a
  category leader for structural homology. Adopt via FFI; only rewrite if
  embedding it as a library becomes important.
- Docking is mostly empirical scoring + search. AutoDock Vina is BSD-3,
  small, and a reasonable Rust port target. GNINA and DiffDock are CNN-based;
  same calculus as the structure-prediction tools.
- License watch: MaxQuant **restricted (free use, not OSI)**, MSFragger
  **academic-only**, DIA-NN **commercial from 1.9.2**, Spectronaut **commercial**,
  Glide **commercial (Schrödinger)**, AlphaFold 2 **Apache-2.0** (weights CC BY-NC),
  AlphaFold 3 **CC BY-NC-SA-4.0**, RoseTTAFold MIT (weights non-commercial),
  most analysis libraries permissive.
