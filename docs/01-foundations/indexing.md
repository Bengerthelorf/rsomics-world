# Indexing

> Random-access indexes for compressed bioinformatics files: fai, bai, csi,
> tbi, gzi.

## Scope

The sidecar index files that turn a streaming BGZF/gzipped archive into a
random-access database: FASTA index (`.fai`), BAM index (`.bai`), coordinate
sort index (`.csi`), tabix index (`.tbi`), and BGZF virtual-offset index
(`.gzi`). Includes the `tabix` CLI that builds them. The codecs themselves
(BGZF / gzip) live in [`compression.md`](compression.md); record-level
readers that *use* these indexes live in [`io-formats.md`](io-formats.md).

## Design notes

- All five formats are well-specified in the
  [hts-specs](http://samtools.github.io/hts-specs/) repository, and all five
  are already implemented in `noodles-*` crates. The Rust gap here is
  near-zero.
- CSI generalises BAI (32-bit → 64-bit coordinates, configurable bin
  depth). Most new pipelines should default to CSI; BAI is kept for
  compatibility with chromosomes ≤ 512 Mbp.
- `.gzi` is a small but essential index for random access into a *gzipped*
  (non-BGZF) FASTA. `samtools faidx` will produce both `.fai` and `.gzi`
  for `.fa.gz`.
- The work is in producing **bit-identical indexes** to htslib so that
  downstream tools (IGV, GATK, bcftools) accept them without complaint. Our
  benchmark plan calls for byte-comparison against `samtools index` output.
- A multi-threaded CSI/BAI builder is a real performance opportunity:
  `samtools index -@` parallelises only BGZF decompression, not the
  bin-tree construction.

## TODO

- [x] **`.fai` FASTA index** — record-offset table for random access into a FASTA file.
  - Reference impl: `C` · [samtools/htslib/faidx.c](https://github.com/samtools/htslib) · `MIT`
  - Existing Rust: [`noodles-fasta`](https://crates.io/crates/noodles-fasta) `0.61.0` (fai reader + writer); `rust-htslib::faidx` (FFI)
  - Existing Rust kind: `pure-port` (noodles) / `FFI-wrapper` (rust-htslib)
  - Existing non-C alternatives: `pyfaidx` (Python)
  - Parallelism: single-threaded (index build is trivially small)
  - SIMD: none needed (linear scan)
  - Quadrant: ① (noodles)
  - GPU-amenable: no — index build is sub-second, no compute upside
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Add CI test that builds `.fai` for GRCh38 and diffs byte-by-byte against `samtools faidx` output.

- [x] **`.bai` BAM index** — BAM coordinate-sort index (≤ 512 Mbp chromosomes).
  - Reference impl: `C` · [samtools/htslib](https://github.com/samtools/htslib) · `MIT`
  - Existing Rust: [`noodles-bam`](https://crates.io/crates/noodles-bam) `0.89.0` (reader + writer)
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: `htsjdk` (Java)
  - Parallelism: serial index build today; a rayon-parallel bin-tree builder is open work
  - SIMD: none
  - Quadrant: ①
  - GPU-amenable: no — irregular tree construction, latency-dominated
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Watch for edge cases in placed-unmapped reads and pseudo-bin 37450 (the BAI metadata bin) — usual cause of incompatibility with older readers.

- [x] **`.csi` Coordinate Sort Index** — 64-bit-capable generalisation of BAI/TBI.
  - Reference impl: `C` · [samtools/htslib](https://github.com/samtools/htslib) · `MIT`
  - Existing Rust: [`noodles-csi`](https://crates.io/crates/noodles-csi) `0.56.0`
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: `htsjdk`
  - Parallelism: serial build; same parallel-builder opportunity as BAI
  - SIMD: none
  - Quadrant: ①
  - GPU-amenable: no — irregular tree construction
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Default to CSI v1 with `min_shift=14` for new outputs. Required for plant genomes (some chromosomes > 512 Mbp).

- [x] **`.tbi` tabix index** — generic positional index for any tab-delimited bgzipped file (BED, GFF, VCF text).
  - Reference impl: `C` · [samtools/htslib](https://github.com/samtools/htslib) · `MIT`
  - Existing Rust: [`noodles-tabix`](https://crates.io/crates/noodles-tabix) `0.62.0`
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: `htsjdk`
  - Parallelism: serial build
  - SIMD: none
  - Quadrant: ①
  - GPU-amenable: no — irregular tree construction
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: The CSI variant (`.tbi` → `.csi`) is preferred for new work; keep `.tbi` writer for backward compat.

- [x] **`.gzi` BGZF virtual-offset index** — sidecar for random access into gzipped FASTA.
  - Reference impl: `C` · [samtools/htslib/bgzf.c](https://github.com/samtools/htslib) · `MIT`
  - Existing Rust: [`noodles-bgzf`](https://crates.io/crates/noodles-bgzf) `0.47.0` (read + write); consumed by `noodles-fasta` for bgzipped references
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: —
  - Parallelism: single-threaded
  - SIMD: none
  - Quadrant: ①
  - GPU-amenable: no — trivial offset table
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Trivial format but a frequent source of confusion. Document the "build `.gzi` *first*, then `.fai`" sequence explicitly in the `rsomics-faidx` CLI. The companion `rsomics-zip` binary ([`compression.md`](compression.md)) needs to learn `.gzi` emission alongside BGZF output — first follow-up TODO over there.

- [ ] **`tabix` (CLI)** — command-line index builder.
  - Reference impl: `C` · [samtools/htslib](https://github.com/samtools/htslib) · `MIT`
  - Existing Rust: no published `tabix`-equivalent CLI; building blocks live in `noodles-tabix` + `noodles-csi`
  - Existing Rust kind: `none`
  - Existing non-C alternatives: —
  - Parallelism: opportunity — rayon-parallel bin-tree construction (samtools serialises this)
  - SIMD: none
  - Quadrant: —
  - GPU-amenable: no — irregular tree construction
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `B` (tool — `rsomics-tabix`)
  - Consumes primitives: `noodles-tabix`, `noodles-csi`, `noodles-bgzf`
  - Notes: Real opportunity. Thin wrapper over noodles primitives but worth shipping for parallel-index-building (tracked upstream as [samtools/htslib#1735](https://github.com/samtools/htslib/issues/1735)). Match htslib output byte-for-byte.
