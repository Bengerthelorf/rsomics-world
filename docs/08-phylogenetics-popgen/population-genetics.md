# Population genetics

> Variant-based population structure, admixture, ancestry inference, IBD,
> low-coverage handling, GWAS support, and scalable genotype data analysis.

## Scope

Includes: classical PLINK-family pipelines (PLINK1.9, PLINK2), ancestry
inference (ADMIXTURE, fastSTRUCTURE, RFMix), PCA / smartPCA (EIGENSOFT),
low-coverage genotype-likelihood frameworks (ANGSD), IBD segment detection
(IBDseq), VCF utilities (vcftools), and scalable array-native popgen
frameworks (sgkit, Hail). Excludes: variant calling itself (see
[../02-genomics](../02-genomics)) and per-sample VCF parsing (covered by
`noodles-vcf`).

## Design notes

- The popgen ecosystem is split between **CLI binaries on flat files**
  (PLINK2 / ADMIXTURE / EIGENSOFT / vcftools) and **array-native modern
  frameworks** (Hail on Spark/Scala, sgkit on Xarray/Dask). The Rust play
  is more interesting in the second camp: `polars` + `arrow` give us a
  natural substrate for sgkit-style population-genetic primitives that
  can be composed in pipelines, without inheriting Python.
- PLINK2's `.pgen` format is heavily-optimized bit-packed genotypes.
  Pure-Rust readers exist but are experimental.
- ANGSD's value is GLs (genotype likelihoods) for low-coverage data;
  algorithmically distinct from PLINK.
- Hail is JVM/Spark — Rust gives us nothing there. sgkit is Python +
  Xarray — Rust can complement by providing fast kernels via PyO3.
- License watch: PLINK2 **GPL-3**, ADMIXTURE **proprietary-free**,
  EIGENSOFT mixed (parts BSD-3, parts proprietary), vcftools **GPL-3**,
  sgkit **Apache-2.0**, Hail **MIT**, ANGSD **GPL-3**, RFMix **GPL-3**,
  IBDseq **proprietary-free**.

## TODO

- [ ] **`PLINK` / `PLINK2`** — GWAS-grade genotype handling and association.
  - Reference impl: `C++` · [chrchang/plink-ng](https://github.com/chrchang/plink-ng) · `GPL-3`
  - Existing Rust: none verified for full PLINK2 functionality; experimental `pgenlib`-style readers may exist
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `Hail` (JVM/Scala), `sgkit` (Python)
  - Parallelism: upstream pthreads
  - SIMD: upstream SSE/AVX
  - Quadrant: —
  - GPU-amenable: maybe — PCA + GLM SIMT-friendly
  - Upstream license: `GPL-3`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-plink`)
  - Consumes primitives: future `rsomics-pgen` (.pgen reader), `noodles-vcf`, `ndarray-linalg`, `polars`, future `rsomics-stats`
  - Notes: PLINK2's value is the `.pgen` data format + speed of common operations (LD pruning, PCA, association tests). Rust port = pgen reader + a focused subcommand set (`--freq`, `--missing`, `--indep-pairwise`, `--glm`). Full feature parity is multi-year.

- [ ] **`ADMIXTURE`** — ML ancestry-proportion inference.
  - Reference impl: `C++` · binary-only from [dalexander.github.io/admixture](https://dalexander.github.io/admixture/) · proprietary-free (no source)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `fastSTRUCTURE` (Python, GPL); `Ohana` (C++)
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — block-relaxation EM is dense
  - Upstream license: proprietary-free (no source)
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-admixture`)
  - Consumes primitives: `ndarray-linalg`, future `rsomics-pgen`, future `rsomics-stats`
  - Notes: Closed source upstream is awkward. The published algorithm (block-relaxation EM on a logistic model) is reproducible from the paper. A clean-room Rust `rsomics-admixture` would be unique in the ecosystem.

- [ ] **`vcftools`** — VCF filtering, statistics, and conversion.
  - Reference impl: `Perl` + `C++` · [vcftools/vcftools](https://github.com/vcftools/vcftools) · `GPL-3`
  - Existing Rust: most common operations already exist via `noodles-vcf` and `vcfexpress`
  - Existing Rust kind: `partial-port` (the common operations have Rust equivalents via noodles + vcfexpress, though no direct port)
  - Existing non-C alternatives: `bcftools` (C, htslib)
  - Parallelism: upstream serial
  - SIMD: none
  - Quadrant: ①
  - GPU-amenable: no — record-level filtering
  - Upstream license: `GPL-3`
  - Priority: `P1`
  - Layer: `subcommand-of-rsomics-bcftools` (cross-listed with [`../02-genomics/variant-calling.md`](../02-genomics/variant-calling.md))
  - Consumes primitives: `noodles-vcf`, `vcfexpress`
  - Notes: Largely superseded by `bcftools` and by `noodles-vcf` + `vcfexpress` (Rust). Adopt the Rust side; no need for a vcftools port per se.

- [ ] **`bcftools`** — cross-reference; canonical entry in [`../02-genomics/variant-calling.md`](../02-genomics/variant-calling.md).
  - Reference impl: `C` · [samtools/bcftools](https://github.com/samtools/bcftools) · `MIT`
  - Existing Rust: see canonical entry
  - Existing Rust kind: see canonical entry
  - Existing non-C alternatives: see canonical entry
  - Parallelism: see canonical entry
  - SIMD: see canonical entry
  - Quadrant: see canonical entry
  - GPU-amenable: see canonical entry
  - Upstream license: `MIT`
  - Priority: cross-listed
  - Layer: `subcommand-of-rsomics-bcftools`
  - Consumes primitives: see canonical entry
  - Notes: Listed here only as a pointer. **Cross-reference only — canonical entry is in `02-genomics/variant-calling.md`.**

- [ ] **`EIGENSOFT`** — smartPCA + tests for ancestry and selection.
  - Reference impl: `C` · [DReichLab/EIG](https://github.com/DReichLab/EIG) · mixed (parts BSD-3, parts proprietary historical components)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: smartPCA-equivalent in `sgkit`, `PLINK2 --pca`
  - Parallelism: upstream pthreads
  - SIMD: BLAS
  - Quadrant: —
  - GPU-amenable: yes — PCA is dense linear algebra
  - Upstream license: mixed (BSD-3 + proprietary parts)
  - Priority: `P2`
  - Layer: `B` (tool — `rsomics-eigensoft`)
  - Consumes primitives: `ndarray-linalg`, future `rsomics-pgen`, future `rsomics-stats`
  - Notes: smartPCA is `nalgebra`-grade Rust math (~few hundred LOC). The valuable extensions (D-statistics, f3/f4, qpAdm) are larger but well-defined. Niche but loved in ancient-DNA / human-popgen communities.

- [ ] **`fastSTRUCTURE`** — variational-Bayes admixture-style inference.
  - Reference impl: `Python` (Cython) · [rajanil/fastStructure](https://github.com/rajanil/fastStructure) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `ADMIXTURE`
  - Parallelism: Python multiprocessing
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — VB updates are dense
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-admixture` (VB-mode alternative to ML mode)
  - Consumes primitives: same as ADMIXTURE plus VB updates in `ndarray-linalg`
  - Notes: Similar product to ADMIXTURE with VB instead of ML. Smaller codebase; could be a thin Rust port.

- [ ] **`sgkit`** — scalable popgen on Xarray + Zarr.
  - Reference impl: `Python` · [sgkit-dev/sgkit](https://github.com/sgkit-dev/sgkit) · `Apache-2.0`
  - Existing Rust: none verified (but `polars` + `zarrs` are the natural Rust counterparts)
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Dask
  - SIMD: BLAS
  - Quadrant: —
  - GPU-amenable: maybe — array-native ops are GPU-friendly via downstream sub-deps
  - Upstream license: `Apache-2.0`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-plink` (sgkit-compatible kernels as a `--backend zarr` mode)
  - Consumes primitives: `zarrs`, `polars`, `ndarray-linalg`
  - Notes: Don't reimplement; instead build Rust kernels that emit / read Zarr stores compatible with sgkit, and expose them as PyO3 wheels. Interop > duplication.

- [ ] **`Hail`** — distributed popgen on Spark/Scala.
  - Reference impl: `Scala` + `Python` · [hail-is/hail](https://github.com/hail-is/hail) · `MIT`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: Spark
  - SIMD: JVM
  - Quadrant: —
  - GPU-amenable: maybe — Spark distribution
  - Upstream license: `MIT`
  - Priority: `P2`
  - Layer: —
  - Consumes primitives: —
  - Notes: JVM/Spark world. No Rust play except possibly via Arrow Flight for hand-off. Skip.

- [ ] **`ANGSD`** — genotype-likelihood framework for low-coverage data.
  - Reference impl: `C++` · [ANGSD/angsd](https://github.com/ANGSD/angsd) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — SAF/SFS estimation is dense
  - Upstream license: `GPL-3`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-angsd`)
  - Consumes primitives: `noodles-bam`, `noodles-vcf`, `ndarray-linalg`, future `rsomics-stats`
  - Notes: Distinct enough algorithmically to justify its own crate (`rsomics-gl` for genotype likelihoods at the foundation layer if shared). The SAF / SFS estimation math is well-bounded.

- [ ] **`RFMix`** — local ancestry inference.
  - Reference impl: `C++` · [slowkoni/rfmix](https://github.com/slowkoni/rfmix) · `GPL-3`
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `Gnomix` (Python neural net)
  - Parallelism: upstream pthreads
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — random-forest scoring
  - Upstream license: `GPL-3`
  - Priority: `P2`
  - Layer: `subcommand-of-rsomics-admixture` (local-ancestry mode)
  - Consumes primitives: `noodles-vcf`, `linfa` (random forest), future `rsomics-haplotype`
  - Notes: Niche but stable. Random-forest based; small ML model. Rust port is feasible with `linfa` once VCF and a haplotype representation are in place.

- [ ] **`IBDseq`** — IBD segment detection on unphased SNPs.
  - Reference impl: `Java` · [Browning lab](https://faculty.washington.edu/sbrowning/ibdseq.html) · proprietary-free (academic-distribution)
  - Existing Rust: none verified
  - Existing Rust kind: `none`
  - Existing non-C alternatives: `hap-ibd` (Java), `iLASH` (C++)
  - Parallelism: JVM threading
  - SIMD: limited
  - Quadrant: —
  - GPU-amenable: maybe — IBD detection over segments parallelises
  - Upstream license: proprietary-free
  - Priority: `P2`
  - Layer: `B` (tool — `rsomics-ibdseq`)
  - Consumes primitives: `noodles-vcf`, `polars`, `rayon`
  - Notes: JVM tool; Rust gives portability and embeddability gains. Algorithm published; small clean-room port possible.
