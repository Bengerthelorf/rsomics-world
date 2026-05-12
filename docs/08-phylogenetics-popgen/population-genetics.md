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
  Pure-Rust readers exist (`pgenlib-rs`, `pgenrust` — verify; not sure
  whether mature). `noodles-vcf` covers VCF/BCF cleanly; we likely want
  a Rust pgen crate to round out the IO layer.
- ANGSD's value is GLs (genotype likelihoods) for low-coverage data;
  algorithmically distinct from PLINK. The C++ codebase is dense but
  well-bounded.
- Hail is JVM/Spark — Rust gives us nothing there. sgkit is Python +
  Xarray — Rust can complement by providing fast kernels via PyO3.
- License watch: PLINK2 **GPL-3**, ADMIXTURE **proprietary-free** (closed
  source, free for academic use), EIGENSOFT mixed (parts BSD-3, parts
  proprietary), vcftools **GPL-3**, sgkit **Apache-2.0**, Hail **MIT**,
  ANGSD **GPL-3**, RFMix **GPL-3** (check), IBDseq **GPL-3**.

## TODO

- [ ] **`PLINK` / `PLINK2`** — GWAS-grade genotype handling and association.
  - Reference impl: `C++` · [chrchang/plink-ng](https://github.com/chrchang/plink-ng) · `GPL-3`
  - Existing Rust: none verified for full PLINK2 functionality; experimental `pgenlib`-style readers may exist
  - Existing non-C alternatives: `Hail` (JVM/Scala), `sgkit` (Python)
  - Priority: `P1`
  - Notes: PLINK2's value is the `.pgen` data format + speed of common
    operations (LD pruning, PCA, association tests). Rust port = pgen
    reader + a focused subcommand set (`rsomics-plink` covering
    `--freq`, `--missing`, `--indep-pairwise`, `--glm`). Full feature parity
    is multi-year.

- [ ] **`ADMIXTURE`** — ML ancestry-proportion inference.
  - Reference impl: `C++` · binary-only from [dalexander.github.io/admixture](https://dalexander.github.io/admixture/) · proprietary-free (no source)
  - Existing Rust: none verified
  - Existing non-C alternatives: `fastSTRUCTURE` (Python, GPL); `Ohana` (C++)
  - Priority: `P1`
  - Notes: Closed source upstream is awkward. The published algorithm
    (block-relaxation EM on a logistic model) is reproducible from the
    paper. A clean-room Rust `rsomics-admixture` would be unique in the
    ecosystem (no other open implementation outside fastSTRUCTURE).

- [ ] **`vcftools`** — VCF filtering, statistics, and conversion.
  - Reference impl: `Perl` + `C++` · [vcftools/vcftools](https://github.com/vcftools/vcftools) · `GPL-3`
  - Existing Rust: most common operations already exist via `noodles-vcf` and `vcfexpress`
  - Existing non-C alternatives: `bcftools` (C, htslib)
  - Priority: `P1`
  - Notes: Largely superseded by `bcftools` (faster, htslib-based) and
    increasingly by `noodles-vcf` + `vcfexpress` (Rust). Adopt the Rust
    side; no need for a vcftools port per se.

- [ ] **`bcftools`** — covered in `02-genomics/variant-calling.md`.
  - Reference impl: `C` · [samtools/bcftools](https://github.com/samtools/bcftools) · `MIT/Expat`
  - Existing Rust: `rust-htslib` (FFI) + `noodles-bcf` (pure-Rust)
  - Existing non-C alternatives: —
  - Priority: cross-listed
  - Notes: Listed here only as a pointer. The variant-calling/IO module
    owns the entry.

- [ ] **`EIGENSOFT`** — smartPCA + tests for ancestry and selection.
  - Reference impl: `C` · [DReichLab/EIG](https://github.com/DReichLab/EIG) · mixed (parts BSD-3, parts proprietary historical components)
  - Existing Rust: none verified
  - Existing non-C alternatives: smartPCA-equivalent in `sgkit`, `PLINK2 --pca`
  - Priority: `P2`
  - Notes: smartPCA is `nalgebra`-grade Rust math (~few hundred LOC). The
    valuable extensions (D-statistics, f3/f4, qpAdm) are larger but well-
    defined. Niche but loved in ancient-DNA / human-popgen communities.

- [ ] **`fastSTRUCTURE`** — variational-Bayes admixture-style inference.
  - Reference impl: `Python` (Cython) · [rajanil/fastStructure](https://github.com/rajanil/fastStructure) · `GPL-3` (per repo)
  - Existing Rust: none verified
  - Existing non-C alternatives: `ADMIXTURE`
  - Priority: `P2`
  - Notes: Similar product to ADMIXTURE with VB instead of ML. Smaller
    codebase; could be a thin Rust port. Skip unless variational-vs-ML
    matters to the user.

- [ ] **`sgkit`** — scalable popgen on Xarray + Zarr.
  - Reference impl: `Python` · [sgkit-dev/sgkit](https://github.com/sgkit-dev/sgkit) · `Apache-2.0`
  - Existing Rust: none verified (but `polars` + `zarr-rs` are the natural Rust counterparts)
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Don't reimplement; instead build Rust kernels that emit / read
    Zarr stores compatible with sgkit, and expose them as PyO3 wheels.
    Interop > duplication.

- [ ] **`Hail`** — distributed popgen on Spark/Scala.
  - Reference impl: `Scala` + `Python` · [hail-is/hail](https://github.com/hail-is/hail) · `MIT`
  - Existing Rust: none
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: JVM/Spark world. No Rust play except possibly via Arrow Flight
    for hand-off. Skip.

- [ ] **`ANGSD`** — genotype-likelihood framework for low-coverage data.
  - Reference impl: `C++` · [ANGSD/angsd](https://github.com/ANGSD/angsd) · `GPL-3`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Distinct enough algorithmically to justify its own crate
    (`rsomics-gl` for genotype likelihoods). The SAF / SFS estimation
    math is well-bounded.

- [ ] **`RFMix`** — local ancestry inference.
  - Reference impl: `C++` · [slowkoni/rfmix](https://github.com/slowkoni/rfmix) · `GPL-3` (check repo)
  - Existing Rust: none verified
  - Existing non-C alternatives: `Gnomix` (Python neural net)
  - Priority: `P2`
  - Notes: Niche but stable. Random-forest based; small ML model. Rust
    port is feasible with `linfa` once VCF and a haplotype representation
    are in place.

- [ ] **`IBDseq`** — IBD segment detection on unphased SNPs.
  - Reference impl: `Java` · [Browning lab](https://faculty.washington.edu/sbrowning/ibdseq.html) · check distribution (proprietary-free)
  - Existing Rust: none verified
  - Existing non-C alternatives: `hap-ibd` (Java), `iLASH` (C++)
  - Priority: `P2`
  - Notes: JVM tool; Rust gives portability and embeddability gains.
    Algorithm published; small clean-room port possible. Lower priority
    than PLINK / ADMIXTURE work.
