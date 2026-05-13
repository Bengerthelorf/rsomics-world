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
  GTDB-Tk **GPL-3**, dRep **MIT**.

## TODO

- [ ] **`MEGAHIT`** — ultra-fast succinct-de-Bruijn-graph metagenome assembler.
  - Reference impl: `C++` · [voutcn/megahit](https://github.com/voutcn/megahit) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads + OpenMP
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — k-mer counting prequel; SdBG construction less so
  - Upstream license: `GPL-3`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-megahit`)
  - Consumes primitives: `rsomics-kmer`, `debruijn` / `ggcat`, `noodles-fastq`, `noodles-fasta`, future succinct dBG crate
  - Notes: The default assembler for short-read metagenomes. Rust port has real value but is non-trivial: succinct de Bruijn graph + multi-k iterative assembly + bubble/tip removal. Probably the largest single rewrite in this module.

- [ ] **`metaSPAdes`** — SPAdes' metagenomic mode.
  - Reference impl: `C++` · [ablab/spades](https://github.com/ablab/spades) · `GPL-2`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE
  - Quadrant: —
  - GPU-amenable: maybe — k-mer counting; graph traversal less so
  - Upstream license: `GPL-2`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-spades` (cross-listed with [`02-genomics/assembly.md`](../02-genomics/assembly.md))
  - Consumes primitives: `rsomics-kmer`, `debruijn`, `noodles-fastq`
  - Notes: Higher contiguity than MEGAHIT, at substantial RAM cost. SPAdes is 100k+ lines of C++ with many sub-pipelines; full port is multi-year. Better to focus on MEGAHIT first and add metaSPAdes-style polishing passes incrementally.

- [ ] **`IDBA-UD`** — iterative de Bruijn assembler for uneven coverage.
  - Reference impl: `C++` · [loneknightpy/idba](https://github.com/loneknightpy/idba) · `GPL-2`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: no — superseded
  - Upstream license: `GPL-2`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Largely superseded by MEGAHIT/metaSPAdes; included for completeness. Skip in favor of MEGAHIT port.

- [ ] **`MetaBAT2`** — adaptive TNF + coverage binner.
  - Reference impl: `C++` · [berkeleylab/metabat](https://bitbucket.org/berkeleylab/metabat) · `BSD-3-Clause-LBNL` (Bitbucket-hosted; gh aliveness N/A)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — TNF computation is SIMT-trivial; label propagation less so
  - Upstream license: `BSD-3-Clause-LBNL`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-metabat`)
  - Consumes primitives: `noodles-fasta`, `noodles-bam`, `ndarray`, `linfa` (clustering), `rayon`
  - Notes: The default short-read binner. Tetranucleotide-frequency distance + abundance graph + label propagation. Self-contained; excellent Rust port target (~few weeks). Bitbucket-hosted upstream.

- [ ] **`MaxBin2`** — EM-based binner with marker genes.
  - Reference impl: `C++` + `Perl` · [sourceforge mirror](https://sourceforge.net/projects/maxbin/) · `BSD`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — EM iteration is dense
  - Upstream license: `BSD`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-metabat` (alternative binning mode)
  - Consumes primitives: same as MetaBAT2 entry
  - Notes: Older method, frequently combined with others under DAS_Tool. Marker-gene EM is the interesting bit; port only after MetaBAT2.

- [ ] **`CONCOCT`** — composition + coverage Gaussian-mixture binner.
  - Reference impl: `Python` (NumPy/scikit-learn) · [BinPro/CONCOCT](https://github.com/BinPro/CONCOCT) · `FreeBSD`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: BLAS
  - Quadrant: —
  - GPU-amenable: maybe — GMM is GPU-friendly
  - Upstream license: `FreeBSD` (BSD-2-Clause)
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-metabat`
  - Consumes primitives: `linfa-clustering` (GMM), `ndarray-linalg`
  - Notes: Python wrapper around scikit-learn GMM. Trivial to reimplement with `linfa`, but value is limited — MetaBAT2/SemiBin2 outperform it.

- [ ] **`SemiBin2`** — self-supervised deep-learning binner.
  - Reference impl: `Python` (PyTorch) · [BigDataBiology/SemiBin](https://github.com/BigDataBiology/SemiBin) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: PyTorch GPU
  - SIMD: PyTorch kernels
  - Quadrant: —
  - GPU-amenable: yes — siamese network inference is dense DL
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-metabat` (DL-binner mode)
  - Consumes primitives: `candle` or `burn`, `noodles-fasta`, `ndarray`
  - Notes: Siamese network on contig features + must-link constraints. Inference with `candle`/`burn` is feasible; training stays PyTorch. Most of the value is the trained model.

- [ ] **`VAMB`** — variational-autoencoder binner.
  - Reference impl: `Python` (PyTorch) · [RasmussenLab/vamb](https://github.com/RasmussenLab/vamb) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: PyTorch GPU
  - SIMD: PyTorch kernels
  - Quadrant: —
  - GPU-amenable: yes — VAE inference is dense DL
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-metabat`
  - Consumes primitives: same as SemiBin2 entry
  - Notes: Same shape as SemiBin2 — VAE/AAE training in PyTorch, but inference and clustering is portable. TaxVamb adds taxonomy semi-supervision; port the inference + clustering only.

- [ ] **`DAS_Tool`** — bin-refinement consensus across multiple binners.
  - Reference impl: `R` + `C++` · [cmks/DAS_Tool](https://github.com/cmks/DAS_Tool) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: R BiocParallel
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: no — greedy bin selection, latency-bound
  - Upstream license: `BSD-3-Clause`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-metabat` (consensus mode)
  - Consumes primitives: future `rsomics-prodigal` (gene caller), `polars`, `linfa`
  - Notes: Self-contained marker-gene-aware bin-refinement greedy algorithm. Easy Rust port if `prodigal`/`pyrodigal`-equivalent is available. The R layer is just glue.

- [ ] **`CheckM`** v1 — marker-gene-based MAG completeness/contamination QC.
  - Reference impl: `Python` · [Ecogenomics/CheckM](https://github.com/Ecogenomics/CheckM) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: no — superseded
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: Largely superseded by CheckM2. Skip in favor of CheckM2 port.

- [ ] **`CheckM2`** — ML-based MAG QC.
  - Reference impl: `Python` (TF/Keras + LightGBM) · [chklovski/CheckM2](https://github.com/chklovski/CheckM2) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: PyTorch / LightGBM threading
  - SIMD: TF/LightGBM kernels
  - Quadrant: —
  - GPU-amenable: yes — GBM inference + protein prediction
  - Upstream license: `GPL-3`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-checkm`)
  - Consumes primitives: future `rsomics-prodigal`, `candle` or LightGBM-Rust binding, ONNX-runtime crate
  - Notes: Default MAG QC tool in 2025. Inference path: Prodigal proteins → KO annotation → feature vector → GBM. Rust port = `pyrodigal`-equiv + `lightgbm-rs` inference. Models distributed as `.pkl` — convert to a Rust-friendly format (ONNX) for deployment.

- [ ] **`GTDB-Tk`** — GTDB taxonomy assignment toolkit.
  - Reference impl: `Python` (wraps HMMER, pplacer, Mash, FastANI) · [Ecogenomics/GTDBTk](https://github.com/Ecogenomics/GTDBTk) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Python multiprocessing + upstream pthreads
  - SIMD: inherits from HMMER / FastANI
  - Quadrant: —
  - GPU-amenable: maybe — HMMER inner loops have GPU variants
  - Upstream license: `GPL-3`
  - Priority: `P0`
  - Layer: `B` (tool — `rsomics-gtdbtk`)
  - Consumes primitives: `sourmash` (mash substitute), `skani` (FastANI substitute), future `rsomics-hmm`, future `rsomics-pplacer` (or EPA-NG via FFI), `polars`
  - Notes: Biggest blocker is `pplacer` (OCaml, ~unmaintained). Path: replace pplacer with EPA-NG (C++, [Pbdas/epa-ng](https://github.com/Pbdas/epa-ng), GPL-3) — already done in some forks. Then orchestration is Rust + `sourmash` (mash substitute) + HMMER-equivalent.

- [ ] **`dRep`** — genome dereplication and ANI clustering.
  - Reference impl: `Python` · [MrOlm/drep](https://github.com/MrOlm/drep) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `skDER` ([raufs/skDER](https://github.com/raufs/skDER), faster C++/Python successor)
  - Parallelism: Python multiprocessing + Mash/nucmer pthreads
  - SIMD: BLAS via upstream tools
  - Quadrant: —
  - GPU-amenable: no — ANI computation, latency-bound
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-drep`)
  - Consumes primitives: `sourmash`, `skani`, `linfa-clustering`
  - Notes: With `sourmash` + `skani` already Rust, dRep becomes mostly a clustering loop. Trivial.

- [x] **`skani`** — fast ANI for metagenomes (Rust-native).
  - Reference impl: `Rust` · [bluenote-1577/skani](https://github.com/bluenote-1577/skani) · `MIT`
  - Existing Rust: [`skani`](https://crates.io/crates/skani) `0.1.1`
  - Existing Rust kind: `rust-native`
  - Existing non-C alternatives: `FastANI` (C++)
  - Parallelism: rayon
  - SIMD: auto-vectorize on sketch comparison
  - Quadrant: ①
  - GPU-amenable: maybe — sketch comparison parallelises
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: The relevant Rust tool for the ANI step used by GTDB-Tk and dRep. Adopt as-is.
