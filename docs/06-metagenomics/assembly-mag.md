# Metagenome assembly and MAG recovery

> Assembly, binning, bin refinement, dereplication, quality control, and
> taxonomic assignment of metagenome-assembled genomes (MAGs).

## Scope

Includes: short-read and hybrid de novo metagenome assemblers (MEGAHIT,
metaSPAdes, IDBA-UD); composition+coverage binners (MetaBAT2, MaxBin2,
CONCOCT); deep-learning binners (SemiBin2, VAMB); bin refinement (DAS_Tool);
genome QC (CheckM, CheckM2); GTDB taxonomy assignment (GTDB-Tk); and genome
dereplication (dRep). Excludes: long-read-only assemblers covered under
genomics (Flye, hifiasm-meta) and read classification (see
[classification](classification.md)).

## Design notes

- Metagenome assembly is a de Bruijn graph problem with heavy memory pressure.
  MEGAHIT's succinct de Bruijn graph (SdBG) is the state of the art for
  RAM/throughput trade-off; metaSPAdes pays more RAM for better contiguity.
  Pure-Rust SdBG is a substantial undertaking but Rust's ownership model
  fits the iterative graph-simplification passes well.
- Binning has bifurcated: classic TNF + coverage clustering (MetaBAT2,
  MaxBin2, CONCOCT) is still useful and easy to port; deep-learning binners
  (SemiBin2, VAMB) require PyTorch-equivalent inference at minimum. `candle`
  and `burn` are realistic for inference; training stays on PyTorch.
- CheckM2 is a gradient-boosted classifier (LightGBM/XGBoost) on top of
  Prodigal protein predictions. The inference path is ~50 lines of NumPy
  in the reference; we need a Rust GBM inference path (`lightgbm-rs`,
  `gbdt-rs` exist).
- GTDB-Tk is essentially `pplacer`/`HMMER`/`mash` wrapped in Python. The
  Rust win is the orchestration + `mash`-equivalent (sourmash works);
  `pplacer` itself is OCaml and is the hardest dependency to displace.
- dRep is a Python pipeline around `Mash`, `nucmer`, and clustering. With
  `sourmash` and a Rust `nucmer` equivalent, a clean-room rewrite is small.
- License watch: MEGAHIT **GPL-3**, SPAdes **GPL-2**, IDBA **GPL-2**,
  MetaBAT2 **BSD-3-Clause-LBNL**, MaxBin2 **BSD-3**, CONCOCT **FreeBSD**,
  SemiBin **MIT**, VAMB **MIT**, DAS_Tool **BSD-3**, CheckM/CheckM2 **GPL-3**,
  GTDB-Tk **GPL-3**, dRep **MIT**. Most assemblers and CheckM are GPL —
  port carefully or rewrite clean-room.

## TODO

- [ ] **`MEGAHIT`** — ultra-fast succinct-de-Bruijn-graph metagenome assembler.
  - Reference impl: `C++` · [voutcn/megahit](https://github.com/voutcn/megahit) · `GPL-3`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: The default assembler for short-read metagenomes. Rust port has
    real value but is non-trivial: succinct de Bruijn graph + multi-k
    iterative assembly + bubble/tip removal. Probably the largest single
    rewrite in this module. Stage as `rsomics-megahit` with the SdBG as a
    standalone reusable crate.

- [ ] **`metaSPAdes`** — SPAdes' metagenomic mode.
  - Reference impl: `C++` · [ablab/spades](https://github.com/ablab/spades) · `GPL-2`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Higher contiguity than MEGAHIT, at substantial RAM cost. SPAdes
    is 100k+ lines of C++ with many sub-pipelines; full port is multi-year.
    Better to focus on MEGAHIT first and add metaSPAdes-style polishing
    passes incrementally.

- [ ] **`IDBA-UD`** — iterative de Bruijn assembler for uneven coverage.
  - Reference impl: `C++` · [loneknightpy/idba](https://github.com/loneknightpy/idba) · `GPL-2`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Largely superseded by MEGAHIT/metaSPAdes; included for
    completeness. Unmaintained upstream. Skip in favor of MEGAHIT port.

- [ ] **`MetaBAT2`** — adaptive TNF + coverage binner.
  - Reference impl: `C++` · [bitbucket berkeleylab/metabat](https://bitbucket.org/berkeleylab/metabat) · `BSD-3-Clause-LBNL`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: The default short-read binner. Tetranucleotide-frequency
    distance + abundance graph + label propagation. Self-contained;
    excellent Rust port target (~few weeks). `ndarray` + `rayon`.

- [ ] **`MaxBin2`** — EM-based binner with marker genes.
  - Reference impl: `C++` + `Perl` · [Linked from bioconda; sourceforge](https://sourceforge.net/projects/maxbin/) · `BSD`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Older method, frequently combined with others under DAS_Tool.
    Marker-gene EM is the interesting bit; port only after MetaBAT2.

- [ ] **`CONCOCT`** — composition + coverage Gaussian-mixture binner.
  - Reference impl: `Python` (NumPy/scikit-learn) · [BinPro/CONCOCT](https://github.com/BinPro/CONCOCT) · `FreeBSD (BSD-2-Clause-like)`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Python wrapper around scikit-learn GMM. Trivial to reimplement
    with `linfa` (Rust ML), but value is limited — MetaBAT2/SemiBin2 outperform it.

- [ ] **`SemiBin2`** — self-supervised deep-learning binner.
  - Reference impl: `Python` (PyTorch) · [BigDataBiology/SemiBin](https://github.com/BigDataBiology/SemiBin) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Siamese network on contig features + must-link constraints.
    Inference with `candle`/`burn` is feasible; training stays PyTorch.
    Most of the value is the trained model; Rust gives us a deployable
    inference binary independent of the Python ecosystem.

- [ ] **`VAMB`** — variational-autoencoder binner.
  - Reference impl: `Python` (PyTorch) · [RasmussenLab/vamb](https://github.com/RasmussenLab/vamb) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Same shape as SemiBin2 — VAE/AAE training in PyTorch, but
    inference and clustering is portable. TaxVamb adds taxonomy
    semi-supervision; port the inference + clustering only.

- [ ] **`DAS_Tool`** — bin-refinement consensus across multiple binners.
  - Reference impl: `R` + `C++` · [cmks/DAS_Tool](https://github.com/cmks/DAS_Tool) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Self-contained marker-gene-aware bin-refinement greedy algorithm.
    Easy Rust port if `prodigal`/`pyrodigal`-equivalent is available. The R
    layer is just glue.

- [ ] **`CheckM`** v1 — marker-gene-based MAG completeness/contamination QC.
  - Reference impl: `Python` · [Ecogenomics/CheckM](https://github.com/Ecogenomics/CheckM) · `GPL-3`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Largely superseded by CheckM2. Skip in favor of CheckM2 port.

- [ ] **`CheckM2`** — ML-based MAG QC (gradient-boosted on Prodigal proteins).
  - Reference impl: `Python` (TF/Keras + LightGBM) · [chklovski/CheckM2](https://github.com/chklovski/CheckM2) · `GPL-3`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Default MAG QC tool in 2025. Inference path: Prodigal proteins
    → KO annotation → feature vector → GBM. Rust port = `pyrodigal`-equiv
    + `lightgbm-rs` inference. Models distributed as `.pkl` — convert to
    a Rust-friendly format (`onnx`?) for deployment.

- [ ] **`GTDB-Tk`** — GTDB taxonomy assignment toolkit.
  - Reference impl: `Python` (wraps HMMER, pplacer, Mash, FastANI) · [Ecogenomics/GTDBTk](https://github.com/Ecogenomics/GTDBTk) · `GPL-3`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Biggest blocker is `pplacer` (OCaml, ~unmaintained). Path:
    replace pplacer with EPA-NG (C++, [Pbdas/epa-ng](https://github.com/Pbdas/epa-ng), GPL-3) — already done in some forks. Then orchestration is Rust + `sourmash` (mash substitute) + HMMER-equivalent. Huge user base; major win.

- [ ] **`dRep`** — genome dereplication and ANI clustering.
  - Reference impl: `Python` · [MrOlm/drep](https://github.com/MrOlm/drep) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: `skDER` (faster C++/Python successor; [raufs/skDER](https://github.com/raufs/skDER))
  - Priority: `P1`
  - Notes: With `sourmash` + a Rust `nucmer`-equivalent (open need —
    consider porting `skani` which is already Rust at [bluenote-1577/skani](https://github.com/bluenote-1577/skani) ★), dRep becomes mostly a clustering loop. Trivial.

- [x] **`skani`** — fast ANI for metagenomes (already Rust).
  - Reference impl: `Rust` · [bluenote-1577/skani](https://github.com/bluenote-1577/skani) · `MIT`
  - Existing Rust: same as reference
  - Existing non-C alternatives: `FastANI` (C++)
  - Priority: `P0` (adopt)
  - Notes: Not in the original task list but it is the relevant Rust tool
    for the ANI step used by GTDB-Tk and dRep. Adopt as-is.
