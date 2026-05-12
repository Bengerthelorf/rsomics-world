# Structure analysis & comparison

> PDB/mmCIF/MD-trajectory IO, structure visualization, secondary-structure
> and accessibility computation, and structural homology / alignment.

## Scope

Includes: PDB/mmCIF/PDBx parsers and structure-data libraries (BioPython.PDB,
MDAnalysis, ProDy), visualization (PyMOL, ChimeraX), descriptors (DSSP,
freesasa), and structural homology / alignment (Foldseek, TM-align, FATCAT).
Excludes: structure prediction (see [structure-prediction](structure-prediction.md))
and docking (see [docking](docking.md)).

## Design notes

- PDB/mmCIF parsing is the equivalent of FASTA/BAM IO for structure work.
  Rust has `pdbtbx` (formerly `rust-pdb`) and a few mmCIF crates; need
  a single canonical pick. `pdbtbx` is the leading candidate.
- Visualization is GUI work (PyMOL, ChimeraX) — not a Rust target. We
  link to them and ship structures that load cleanly.
- Foldseek is the category-leading structural-homology tool and is
  already hand-optimized C++/SIMD. We adopt via FFI (or call as a
  process) rather than rewrite.
- TM-align and FATCAT are small, focused, classical alignment kernels.
  TM-align in particular is ~3k lines of C++ and is a clean port target.
- DSSP and freesasa are well-bounded geometric kernels (~1-3k LOC each).
  Both worth pure-Rust ports — small surface area, big breadth-of-use.
- MD-trajectory analysis (MDAnalysis) is a much larger undertaking — DCD,
  XTC, TRR, AMBER, NAMD formats; periodic boundary conditions; selection
  language. A pure-Rust MDAnalysis is a multi-year effort and probably
  out of scope for now.
- License watch: PyMOL **BSD-like open + paid Schrödinger commercial**,
  ChimeraX **free for academic, commercial license needed**, MDAnalysis
  **GPL-2+**, BioPython **Biopython License**, ProDy **MIT**, freesasa
  **MIT**, DSSP **BSD-2-Clause** (mkdssp), Foldseek **GPL-3**, TM-align
  research-license (free for academic use, not OSI), FATCAT **academic**.

## TODO

- [ ] **`PyMOL`** — molecular graphics + scripting (industry standard).
  - Reference impl: `C++` + `Python` · [schrodinger/pymol-open-source](https://github.com/schrodinger/pymol-open-source) · BSD-like open + Schrödinger commercial extensions
  - Existing Rust: none verified
  - Existing non-C alternatives: `ChimeraX`, web-based viewers (Mol*)
  - Priority: `P2`
  - Notes: GUI; out of scope. Link to existing tools.

- [ ] **`ChimeraX`** — modern molecular visualization (successor to UCSF Chimera).
  - Reference impl: `C++` + `Python` · [RBVI/ChimeraX](https://github.com/RBVI/ChimeraX) · free academic, commercial requires license
  - Existing Rust: none verified
  - Existing non-C alternatives: `PyMOL`, `Mol*`
  - Priority: `P2`
  - Notes: GUI; out of scope.

- [ ] **`MDAnalysis`** — Python library for MD trajectory analysis.
  - Reference impl: `Python` + `C` (XDR libs) · [MDAnalysis/mdanalysis](https://github.com/MDAnalysis/mdanalysis) · `GPL-2+`
  - Existing Rust: partial — [`chemfiles-rs`](https://github.com/chemfiles/chemfiles.rs) (FFI to C++ `chemfiles`) provides multi-format trajectory IO
  - Existing non-C alternatives: `chemfiles` (C++, MIT)
  - Priority: `P2`
  - Notes: Large surface area (formats + selection language + analysis).
    Realistic Rust delivery: focused crates on top of `chemfiles` for
    common analyses, not a full MDAnalysis port. GPL-2 inheritance is a
    constraint if porting directly.

- [ ] **`BioPython.PDB`** — PDB/mmCIF parser and structure manipulation.
  - Reference impl: `Python` · [biopython/biopython](https://github.com/biopython/biopython) · Biopython License (BSD-style)
  - Existing Rust: [`pdbtbx`](https://github.com/douweschulte/pdbtbx) (pure Rust PDB/mmCIF reader/writer, MIT)
  - Existing non-C alternatives: —
  - Priority: `P0` (adopt `pdbtbx`)
  - Notes: `pdbtbx` is the Rust answer. Adopt it as the structure IO
    layer for `rsomics-structure`.

- [x] **`pdbtbx`** — pure-Rust PDB / mmCIF IO (foundation).
  - Reference impl: `Rust` · [douweschulte/pdbtbx](https://github.com/douweschulte/pdbtbx) · `MIT`
  - Existing Rust: same as reference
  - Existing non-C alternatives: —
  - Priority: `P0` (adopt)
  - Notes: Foundational. Required for any downstream structure work.

- [ ] **`ProDy`** — protein dynamics analysis (NMA, GNM, ANM).
  - Reference impl: `Python` (NumPy) · [prody/ProDy](https://github.com/prody/ProDy) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Normal-mode analysis is `ndarray` + `nalgebra` eigensolves — a
    self-contained Rust crate is straightforward. Useful but niche.

- [ ] **`freesasa`** — solvent-accessible surface area.
  - Reference impl: `C` · [mittinatten/freesasa](https://github.com/mittinatten/freesasa) · `MIT`
  - Existing Rust: none verified (FFI possible via cc-rs)
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Small, clean, geometric kernel (~3k LOC). Excellent pure-Rust
    port target. Ships well as a tiny crate `rsomics-sasa`.

- [ ] **`DSSP`** — secondary structure assignment from coordinates.
  - Reference impl: `C++` · [PDB-REDO/dssp](https://github.com/PDB-REDO/dssp) · `BSD-2-Clause`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Classic, well-specified algorithm. Pure-Rust port is small and
    high-value (every downstream tool calls DSSP-equivalent at some point).
    `rsomics-dssp` is a candidate "starter" crate.

- [x] **`Foldseek`** — fast structural homology (3Di alphabet + SIMD).
  - Reference impl: `C++` (heavily SIMD-optimized) · [steineggerlab/foldseek](https://github.com/steineggerlab/foldseek) · `GPL-3`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P0` (adopt via FFI/process)
  - Notes: Category leader. The 3Di alphabet idea is published; the SIMD
    kernels are not trivially portable. Use Foldseek as the substrate;
    don't rewrite. GPL-3 means we wrap by subprocess, not link.

- [ ] **`TM-align`** — pairwise structural alignment with TM-score.
  - Reference impl: `C++` (and Fortran) · [Zhang lab](https://zhanggroup.org/TM-align/) · academic-use license (not OSI)
  - Existing Rust: none verified
  - Existing non-C alternatives: `Foldseek` includes TM-align as one of its alignment modes
  - Priority: `P1`
  - Notes: ~3k LOC. Algorithm published in detail. Clean-room pure-Rust
    rewrite is realistic and useful as a permissively-licensed alternative
    (TM-align upstream is not OSI). Pair with `pdbtbx`.

- [ ] **`FATCAT`** — flexible structural alignment.
  - Reference impl: `Java` · [Godzik-Lab/FATCAT-public](https://github.com/godziklab/FATCAT-public) · academic
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Niche. Skip unless flexible alignment is core to a downstream
    project.
