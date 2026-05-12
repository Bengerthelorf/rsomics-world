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
- CSI generalises BAI (32-bit `→` 64-bit coordinates, configurable bin
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

- [x] **`.fai` FASTA index** — record-offset table for random access into a
  FASTA file.
  - Reference impl: `C` · [samtools/htslib/faidx.c](https://github.com/samtools/htslib) · `MIT/Expat`
  - Existing Rust: [`noodles-fasta`](https://crates.io/crates/noodles-fasta)
    (`fai` reader + writer); `rust-htslib::faidx` (FFI)
  - Existing non-C alternatives: `pyfaidx` (Python)
  - Priority: `P0`
  - Notes: Adopt noodles. Add CI test that builds `.fai` for GRCh38 and
    diffs against `samtools faidx` output.

- [x] **`.bai` BAM index** — BAM coordinate-sort index (≤ 512 Mbp
  chromosomes).
  - Reference impl: `C` · [samtools/htslib](https://github.com/samtools/htslib) · `MIT/Expat`
  - Existing Rust: [`noodles-bam`](https://crates.io/crates/noodles-bam)
    (reader + writer)
  - Existing non-C alternatives: `htsjdk` (Java)
  - Priority: `P0`
  - Notes: Adopt. Watch for edge cases in placed-unmapped reads and
    pseudo-bin 37450 (the BAI metadata bin) — these are the usual cause of
    incompatibility with older readers.

- [x] **`.csi` Coordinate Sort Index** — 64-bit-capable generalisation of
  BAI/TBI.
  - Reference impl: `C` · [samtools/htslib](https://github.com/samtools/htslib) · `MIT/Expat`
  - Existing Rust: [`noodles-csi`](https://crates.io/crates/noodles-csi)
  - Existing non-C alternatives: `htsjdk`
  - Priority: `P0`
  - Notes: Adopt and default to CSI v1 with `min_shift=14` for new outputs.
    Required for plant genomes (some chromosomes > 512 Mbp).

- [x] **`.tbi` tabix index** — generic positional index for any tab-
  delimited bgzipped file (BED, GFF, VCF text).
  - Reference impl: `C` · [samtools/htslib](https://github.com/samtools/htslib) · `MIT/Expat`
  - Existing Rust: [`noodles-tabix`](https://crates.io/crates/noodles-tabix)
  - Existing non-C alternatives: `htsjdk`
  - Priority: `P0`
  - Notes: Adopt. The CSI variant (`.tbi` → `.csi`) is preferred for new
    work; keep `.tbi` writer for backward compat.

- [x] **`.gzi` BGZF virtual-offset index** — sidecar for random access into
  gzipped FASTA.
  - Reference impl: `C` · [samtools/htslib/bgzf.c](https://github.com/samtools/htslib) · `MIT/Expat`
  - Existing Rust: [`noodles-bgzf`](https://crates.io/crates/noodles-bgzf)
    (read + write); used by `noodles-fasta` for bgzipped references
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Adopt. Trivial format but a frequent source of confusion.
    Document the "build `.gzi` *first*, then `.fai`" sequence explicitly.

- [~] **`tabix` (CLI)** — command-line index builder.
  - Reference impl: `C` · [samtools/htslib](https://github.com/samtools/htslib) · `MIT/Expat`
  - Existing Rust: no published `tabix`-equivalent CLI yet; building blocks
    live in `noodles-tabix` + `noodles-csi`
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Real opportunity. A `rsomics-tabix` CLI on top of
    `noodles-{tabix,csi,bgzf}` is a thin wrapper but worth shipping for
    parallel-index-building (tracked under
    [samtools/htslib#1735](https://github.com/samtools/htslib/issues/1735)
    upstream). Match htslib output byte-for-byte.
