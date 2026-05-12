# Data structures

> The string-index, hashing, and probabilistic structures underlying every
> aligner, k-mer counter, and sketcher in the rsomics stack.

## Scope

Sequence indexes (FM-index, BWT, suffix array), k-mer hashing (ntHash,
MurmurHash3, xxHash), MinHash sketching, HyperLogLog cardinality
estimation, and Bloom / Cuckoo membership filters. Sequence *alignment*
algorithms that consume these indexes live in module 02; *codec-level*
compression lives in [`compression.md`](compression.md).

## Design notes

- The classical full-text indexes (FM-index, BWT, suffix array) are well
  covered by [`rust-bio`](https://github.com/rust-bio/rust-bio) — but the
  implementation is dated and several issues (terminal-sentinel
  requirements, lack of FMD index variants, limited 64-bit support) need
  resolution before we declare it production-ready for a BWA-class
  aligner.
- Hash function choice is performance-critical for any k-mer-based tool.
  `ntHash` (rolling) is the right default for adjacent k-mers; `xxHash3`
  for general use; `MurmurHash3` for compatibility with existing sketches
  (Mash, sourmash, finch).
- MinHash and HyperLogLog have multiple Rust implementations with
  overlapping but inconsistent APIs. `sourmash` (Rust core) and `finch`
  are the production-grade options; new work should plug into one of them
  rather than start a third sketching library.
- For Bloom / Cuckoo filters,
  [`probabilistic-collections`](https://github.com/jeffrey-xiao/probabilistic-collections-rs)
  is the most complete single crate; lots of one-off implementations exist
  but lose to it on either API or performance.
- SIMD matters everywhere here: the inner loops of FM-index rank queries,
  rolling hashes, and Bloom filter lookups are all vectorisable, but only
  a few existing crates expose explicit SIMD paths. Opportunity for a
  `rsomics-kmer` crate that consolidates the fast paths.

## TODO

- [~] **FM-index** — succinct full-text index for backward search.
  - Reference impl: `C++` · [Ferragina & Manzini original; bwa internal](https://github.com/lh3/bwa) · various
  - Existing Rust: [`rust-bio::data_structures::fmindex`](https://docs.rs/bio/latest/bio/data_structures/fmindex/index.html);
    [`fm-index`](https://docs.rs/fm-index/latest/fm_index/);
    [`nucleic-acid`](https://lib.rs/crates/nucleic-acid)
  - Existing non-C alternatives: `sdsl-lite` (C++)
  - Priority: `P0`
  - Notes: rust-bio implementation works but issues
    [#30](https://github.com/rust-bio/rust-bio/issues/30) and
    [#495](https://github.com/rust-bio/rust-bio/issues/495) flag UX +
    correctness corners. A modernised `rsomics-fmindex` with SIMD rank,
    64-bit suffix arrays, and FMD support (needed for BWA-style aligners)
    is on the critical path.

- [~] **BWT (Burrows-Wheeler Transform)** — string permutation underlying
  FM-index.
  - Reference impl: `C` · [bwa BWT routines](https://github.com/lh3/bwa) · `MIT`
  - Existing Rust: `rust-bio::data_structures::bwt`; `nucleic-acid`
  - Existing non-C alternatives: `libdivsufsort` (C/C++)
  - Priority: `P0`
  - Notes: Construction performance is the bottleneck for indexing
    multi-Gbp references. The induced-sorting (SA-IS) variant is the
    state of the art; need a pure-Rust port that matches `libdivsufsort`
    on GRCh38.

- [~] **Suffix array** — sorted suffix offsets.
  - Reference impl: `C` · [libdivsufsort](https://github.com/y-256/libdivsufsort) · `MIT`
  - Existing Rust: [`rust-bio::data_structures::suffix_array`](https://github.com/rust-bio/rust-bio/blob/master/src/data_structures/suffix_array.rs);
    [`suffix`](https://crates.io/crates/suffix);
    [`divsufsort`](https://crates.io/crates/divsufsort) (Rust port of libdivsufsort)
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: `divsufsort` is the right adoption target for SA construction.
    Validate against `libdivsufsort` on real genomes.

- [x] **`ntHash`** — rolling hash for DNA k-mers.
  - Reference impl: `C++` · [bcgsc/ntHash](https://github.com/bcgsc/ntHash) · `MIT`
  - Existing Rust: [`nthash`](https://crates.io/crates/nthash);
    [`nthash_rs`](https://docs.rs/nthash-rs/latest/nthash_rs/) (idiomatic pure-Rust port)
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Adopt. `nthash_rs` is the cleaner modern port (handles
    non-ACGT bases, canonical k-mers). Used by sourmash, GGCAT, many
    others.

- [x] **`MurmurHash3`** — general non-cryptographic hash.
  - Reference impl: `C++` · [aappleby/smhasher](https://github.com/aappleby/smhasher) · `Public domain`
  - Existing Rust: [`murmurhash3`](https://crates.io/crates/murmurhash3);
    [`mur3`](https://docs.rs/mur3);
    [`murmur3`](https://github.com/stusmall/murmur3)
  - Existing non-C alternatives: `xxHash` (faster modern alternative)
  - Priority: `P1`
  - Notes: Adopt one (`mur3` has the cleanest Hasher API). Required for
    compatibility with Mash/sourmash sketches; new internal hashes should
    prefer `xxh3` or `ahash`.

- [x] **MinHash sketching** — locality-sensitive sketch for similarity.
  - Reference impl: `C++` · [marbl/Mash](https://github.com/marbl/Mash) · `BSD-3-Clause`
  - Existing Rust: [`sourmash`](https://crates.io/crates/sourmash) (Rust core),
    [`finch`](https://github.com/onecodex/finch-rs)
  - Existing non-C alternatives: `mash` itself (C++)
  - Priority: `P0`
  - Notes: Adopt sourmash (broader feature set, scaled MinHash) or finch
    (lighter, faster, no Python interop). Both are production-grade. New
    work plugs in as a feature on these crates.

- [x] **HyperLogLog** — approximate cardinality counter.
  - Reference impl: `C++` · [original Flajolet et al.](https://research.neustar.biz/2012/10/25/sketch-of-the-day-hyperloglog-cornerstone-of-a-big-data-pipeline/) · academic
  - Existing Rust: [`probabilistic-collections`](https://crates.io/crates/probabilistic-collections);
    [`hyperloglog`](https://crates.io/crates/hyperloglog);
    [`amadeus-streaming`](https://crates.io/crates/amadeus-streaming)
  - Existing non-C alternatives: HLL ships in Redis, ClickHouse, etc.
  - Priority: `P1`
  - Notes: Adopt `probabilistic-collections` for one-stop import. SIMD
    HLL is possible but rarely the bottleneck — only worth doing if
    profiling justifies it.

- [x] **Bloom filter** — approximate-membership probabilistic set.
  - Reference impl: `C++` · academic; many implementations · public domain
  - Existing Rust: [`probabilistic-collections`](https://crates.io/crates/probabilistic-collections);
    [`bloom-filters`](https://crates.io/crates/bloom-filters);
    [`fastbloom`](https://crates.io/crates/fastbloom) (SIMD)
  - Existing non-C alternatives: —
  - Priority: `P0`
  - Notes: Adopt `fastbloom` for hot paths (SIMD-accelerated), fall back
    to `probabilistic-collections` for the feature-rich API. Used by
    ABySS, sourmash, many metagenomics tools.

- [x] **Cuckoo filter** — Bloom alternative with deletion + better locality.
  - Reference impl: `C++` · [efficient/cuckoofilter](https://github.com/efficient/cuckoofilter) · `Apache-2.0`
  - Existing Rust: [`cuckoofilter`](https://docs.rs/cuckoofilter);
    [`probabilistic-collections`](https://crates.io/crates/probabilistic-collections);
    [`autoscale_cuckoo_filter`](https://crates.io/crates/autoscale_cuckoo_filter)
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: Adopt. Strong fit for k-mer deduplication where deletion is
    needed (streaming metagenomics).

- [ ] **Compacted de Bruijn graph** — adjacency structure for assembly +
  pangenome.
  - Reference impl: `C++` · [GATB-bcalm](https://github.com/GATB/bcalm) · `MIT`
  - Existing Rust: [`rust-debruijn`](https://github.com/10XGenomics/rust-debruijn);
    [`ggcat`](https://github.com/algbio/ggcat) (compacted + coloured);
    [`rust-mdbg`](https://github.com/ekimb/rust-mdbg) (minimizer-space)
  - Existing non-C alternatives: —
  - Priority: `P1`
  - Notes: `ggcat` is production-grade and the right consolidation target.
    `rust-debruijn` is older but used by 10x internals. Cross-references
    [`02-genomics/assembly.md`](../02-genomics/assembly.md).
