# Variant calling

> Small-variant (SNV + indel) callers — germline, somatic, trio, long-read.

## Scope

Callers for SNVs and small indels: GATK HaplotypeCaller, DeepVariant,
FreeBayes, bcftools call, Strelka2, Octopus, Clair3,
PEPPER-Margin-DeepVariant, DeepTrio. Structural-variant callers live in
[`sv-calling.md`](sv-calling.md). Variant annotation (consequence
prediction, population frequency, pathogenicity) lives in
[`annotation.md`](annotation.md).

## Design notes

- Two algorithmic families dominate: **local haplotype assembly + HMM**
  (GATK HC, Strelka2, Octopus) and **deep-learning pileup classification**
  (DeepVariant, Clair3, PEPPER-DeepVariant, DeepTrio). The HMM family is
  natural Rust territory; the DL family needs `candle`/`burn` for
  inference (we ship models, we do not train them).
- `varlociraptor` is the most credible production Rust caller today —
  uncertainty-aware, latent-variable model, handles small + structural
  variants, FDR-controlled filtering. It is not a drop-in HaplotypeCaller
  replacement but is the reference Rust example.
- bcftools `mpileup | call` remains the lightweight default for many
  pipelines and is *not* equivalent to DeepVariant or HaplotypeCaller —
  it is a Bayesian site-by-site caller, no local assembly. Worth a Rust
  port because the algorithm is simple and the existing C is single-
  threaded.
- For long-read calling, Clair3 (HMM + DL hybrid) and
  PEPPER-Margin-DeepVariant (RNN + DV) are the production standards.
  Both rely on neural network inference; `candle` can run the ONNX-exported
  models.
- GATK HaplotypeCaller is the hardest target: 200k+ LOC of Java, decade
  of accumulated heuristics, integration with VQSR / CNN-filtering. A
  Rust port that matches its sensitivity is a multi-year project.
  Practical near-term goal: a `rsomics-haplotypecaller` that matches
  on Genome-in-a-Bottle truth sets to within 1% F1.

## TODO

- [ ] **GATK `HaplotypeCaller`** — local-assembly + HMM germline caller.
  - Reference impl: `Java` · [broadinstitute/gatk](https://github.com/broadinstitute/gatk) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing non-C alternatives: `DRAGEN-GATK` (Illumina, partly FPGA)
  - Priority: `P0`
  - Notes: Highest-impact, hardest target. Algorithm: assemble local
    haplotypes via dBG, realign reads with PairHMM, genotype via
    Bayesian likelihoods. Each stage maps cleanly to Rust modules;
    PairHMM is SIMD-friendly (see PairHMM optimisations in DRAGEN /
    SmithWatermanIntel).

- [ ] **`DeepVariant`** — CNN pileup classifier.
  - Reference impl: `Python / TensorFlow / C++` · [google/deepvariant](https://github.com/google/deepvariant) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Inference workload. Port models to `candle`/`burn` via ONNX;
    pileup-image generation maps naturally to `noodles-bam` + `ndarray`.
    Training stays in TensorFlow upstream — we just ship inference.

- [ ] **`FreeBayes`** — Bayesian haplotype-based caller.
  - Reference impl: `C++` · [freebayes/freebayes](https://github.com/freebayes/freebayes) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Used heavily in non-human variant calling (plants, pathogens).
    MIT licence, modest codebase (~30k LOC), good Rust target. Algorithm
    is well-described in the original paper.

- [ ] **`bcftools call`** — site-by-site Bayesian caller (germline).
  - Reference impl: `C` · [samtools/bcftools](https://github.com/samtools/bcftools) · `MIT/Expat`
  - Existing Rust: noodles-vcf for IO; no caller logic verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: The simplest caller in this list. `bcftools mpileup | call` is
    still the default for non-human, non-clinical variant calling. A
    `rsomics-bcftools` should match it bit-for-bit on the call subcommand.

- [ ] **`Strelka2`** — germline + somatic small-variant caller.
  - Reference impl: `C++` · [Illumina/strelka](https://github.com/Illumina/strelka) · `GPL-3.0`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Top-tier on tumor-normal somatic accuracy. GPL-3.0 forces
    clean-room re-implementation. Active competitor to Mutect2 (GATK)
    for somatic; algorithm differs (graph-based candidate haplotypes).

- [ ] **`Octopus`** — flexible Bayesian caller (germline / somatic / trio /
  polyploid).
  - Reference impl: `C++` · [luntergroup/octopus](https://github.com/luntergroup/octopus) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Most flexible caller of the bunch (one tool, many modes via
    different priors). MIT licence. Active maintenance. Strong candidate
    for a Rust rewrite — algorithm is well-modularised in the C++.

- [ ] **`Clair3`** — DL + HMM long-read caller.
  - Reference impl: `Python / C++` · [HKU-BAL/Clair3](https://github.com/HKU-BAL/Clair3) · `BSD-3-Clause`
  - Existing Rust: none verified (`Clair3-RNA` extension exists)
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Production long-read caller (ONT + PacBio). Inference via
    `candle` is feasible; the tensor pipeline plus pileup featurisation
    is the bulk of the work.

- [ ] **`PEPPER-Margin-DeepVariant`** — nanopore variant + phasing pipeline.
  - Reference impl: `Python / C++` · [kishwarshafin/pepper](https://github.com/kishwarshafin/pepper) · `MIT`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Three-stage pipeline (PEPPER haplotype-aware reads → Margin
    phasing → DeepVariant). Same `candle` story as DeepVariant +
    Clair3. Margin is C++ and self-contained.

- [ ] **`DeepTrio`** — trio variant caller built on DeepVariant.
  - Reference impl: `Python / TensorFlow` · [google/deepvariant](https://github.com/google/deepvariant) (deeptrio subdir) · `BSD-3-Clause`
  - Existing Rust: none verified
  - Existing non-C alternatives: —
  - Priority: `P2`
  - Notes: Re-uses DeepVariant infrastructure with a trio-aware model.
    Trivial to add once DeepVariant inference exists in Rust.

- [x] **`varlociraptor`** — uncertainty-aware variant caller (Rust).
  - Reference impl: `Rust` · [varlociraptor/varlociraptor](https://github.com/varlociraptor/varlociraptor) · `MIT`
  - Existing Rust: [`varlociraptor`](https://crates.io/crates/varlociraptor)
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Adopt. Production-grade Rust caller for the
    uncertainty-aware niche (heterogeneous samples, tumor scenarios,
    FDR-controlled filtering). Composes well with other callers as a
    post-processor. Already used in `nf-core/sarek`.
