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
  - Existing Rust: [`rust-bio::data_structures::fmindex`](https://crates.io/crates/bio) (in `bio` `3.0.0`); [`fm-index`](https://crates.io/crates/fm-index) `0.3.0`; [`nucleic-acid`](https://crates.io/crates/nucleic-acid) `0.1.1`
  - Existing Rust kind: `partial-port` (no FMD-index variant; rank-query SIMD missing)
  - Existing non-C alternatives: `sdsl-lite` (C++)
  - Parallelism: single-threaded construction and query in all three crates today
  - SIMD: none explicit (rank loops are auto-vectorize candidates)
  - Quadrant: ③
  - GPU-amenable: maybe — rank queries on a fixed index port to GPU but the engineering cost is high
  - Upstream license: various (BWA itself MIT)
  - Priority: `P0`
  - Layer: `A` (foundation — `rsomics-fm-index`)
  - Consumes primitives: —
  - Notes: rust-bio implementation works but issues [#30](https://github.com/rust-bio/rust-bio/issues/30) and [#495](https://github.com/rust-bio/rust-bio/issues/495) flag UX + correctness corners. A modernised `rsomics-fm-index` with `std::simd` rank queries, 64-bit suffix arrays, and FMD support (needed for BWA-style aligners) is on the critical path.

- [~] **BWT (Burrows-Wheeler Transform)** — string permutation underlying FM-index.
  - Reference impl: `C` · [bwa BWT routines](https://github.com/lh3/bwa) · `MIT`
  - Existing Rust: `bio::data_structures::bwt` (in `bio` `3.0.0`); `nucleic-acid`
  - Existing Rust kind: `partial-port`
  - Existing non-C alternatives: `libdivsufsort` (C/C++)
  - Parallelism: single-threaded SA construction (the bottleneck)
  - SIMD: none explicit
  - Quadrant: ③
  - GPU-amenable: maybe — SA-IS variants have GPU implementations in the literature; engineering cost is non-trivial
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `A` (foundation — same crate as FM-index)
  - Consumes primitives: —
  - Notes: Construction performance is the bottleneck for indexing multi-Gbp references. The induced-sorting (SA-IS) variant is the state of the art; need a pure-Rust port that matches `libdivsufsort` on GRCh38. Parallel SA-IS (parlay-style) is a known but unexplored Rust opportunity.

- [~] **Suffix array** — sorted suffix offsets.
  - Reference impl: `C` · [y-256/libdivsufsort](https://github.com/y-256/libdivsufsort) · `MIT`
  - Existing Rust: `bio::data_structures::suffix_array` (in `bio` `3.0.0`); [`suffix`](https://crates.io/crates/suffix) `1.3.0`; [`divsufsort`](https://crates.io/crates/divsufsort) `2.0.0` (Rust port of libdivsufsort)
  - Existing Rust kind: `pure-port` (`divsufsort` is a faithful Rust port)
  - Existing non-C alternatives: —
  - Parallelism: single-threaded
  - SIMD: none explicit
  - Quadrant: ③
  - GPU-amenable: maybe — parallel SA construction is researched but engineering-heavy
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `A` (foundation — same crate as FM-index)
  - Consumes primitives: —
  - Notes: `divsufsort` is the right adoption target for SA construction. Validate against `libdivsufsort` on real genomes byte-for-byte.

- [x] **`ntHash`** — rolling hash for DNA k-mers.
  - Reference impl: `C++` · [bcgsc/ntHash](https://github.com/bcgsc/ntHash) · `MIT`
  - Existing Rust: [`nthash`](https://crates.io/crates/nthash) `0.5.1` (luizirber); [`nthash-rs`](https://crates.io/crates/nthash-rs) `0.1.3` (pure-Rust port)
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: —
  - Parallelism: per-k-mer rolling, trivially parallel over input chunks
  - SIMD: none explicit yet (a candidate for `std::simd` rolling-hash batching)
  - Quadrant: ①
  - GPU-amenable: yes — rolling hash over a large sequence batch is SIMT-friendly (per-position parallelism)
  - Upstream license: `MIT`
  - Priority: `P0`
  - Layer: `A` (foundation — `rsomics-kmer`)
  - Consumes primitives: —
  - Notes: `nthash-rs` is the cleaner modern port (handles non-ACGT bases, canonical k-mers). Used by sourmash, GGCAT, many others.

- [x] **`MurmurHash3`** — general non-cryptographic hash.
  - Reference impl: `C++` · [aappleby/smhasher](https://github.com/aappleby/smhasher) · `Public domain`
  - Existing Rust: [`murmurhash3`](https://crates.io/crates/murmurhash3) `0.0.5`; [`mur3`](https://crates.io/crates/mur3) `0.1.0`; [`murmur3`](https://crates.io/crates/murmur3)
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: `xxHash` (faster modern alternative)
  - Parallelism: stateless per-key
  - SIMD: none explicit
  - Quadrant: ④
  - GPU-amenable: yes — per-key stateless hash is SIMT-trivial (but rarely the bottleneck)
  - Upstream license: `Public domain`
  - Priority: `P1`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Pick one (`mur3` has the cleanest Hasher API). Required for compatibility with Mash/sourmash sketches; new internal hashes should prefer `xxh3` or `ahash`.

- [x] **MinHash sketching** — locality-sensitive sketch for similarity.
  - Reference impl: `C++` · [marbl/Mash](https://github.com/marbl/Mash) · `BSD-3-Clause`
  - Existing Rust: [`sourmash`](https://crates.io/crates/sourmash) `0.22.0` (Rust core); [`finch`](https://github.com/onecodex/finch-rs) `0.6.2`
  - Existing Rust kind: `rust-native` (independent Rust impls of an academic algorithm; Mash is one of several MinHash impls, not a C/C++ upstream being ported)
  - Existing non-C alternatives: `mash` itself (C++)
  - Parallelism: rayon-amenable per-sketch
  - SIMD: none explicit in sketching loops; underlying hashes auto-vectorize
  - Quadrant: ①
  - GPU-amenable: maybe — sketch construction is parallel but I/O bound on the FASTA read
  - Upstream license: `BSD-3-Clause`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Adopt sourmash (broader feature set, scaled MinHash) or finch (lighter, faster, no Python interop). Both are production-grade. New work plugs in as a feature on these crates.

- [x] **HyperLogLog** — approximate cardinality counter.
  - Reference impl: `C++` · [original Flajolet et al.](https://research.neustar.biz/2012/10/25/sketch-of-the-day-hyperloglog-cornerstone-of-a-big-data-pipeline/) · academic
  - Existing Rust: [`probabilistic-collections`](https://crates.io/crates/probabilistic-collections) `0.7.0`; [`hyperloglog`](https://crates.io/crates/hyperloglog) `1.0.3`; [`amadeus-streaming`](https://crates.io/crates/amadeus-streaming) `0.4.3` (SIMD-accelerated)
  - Existing Rust kind: `rust-native` (academic algorithm; multiple independent Rust impls, no canonical C upstream being ported)
  - Existing non-C alternatives: HLL ships in Redis, ClickHouse, etc.
  - Parallelism: per-register merge; rayon-able
  - SIMD: explicit in `amadeus-streaming` (other crates rely on auto-vectorize)
  - Quadrant: ① (`amadeus-streaming`) / ④ (others)
  - GPU-amenable: maybe — only worth doing for very large estimation streams
  - Upstream license: academic / public domain in spirit
  - Priority: `P1`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Adopt `probabilistic-collections` for one-stop import unless profiling justifies the SIMD path.

- [x] **Bloom filter** — approximate-membership probabilistic set.
  - Reference impl: `C++` · academic; many implementations · public domain
  - Existing Rust: [`probabilistic-collections`](https://crates.io/crates/probabilistic-collections) `0.7.0`; [`bloom-filters`](https://crates.io/crates/bloom-filters) `0.1.2`; [`fastbloom`](https://github.com/tomtomwombat/fastbloom) `0.17.0` (`no_std`, concurrent)
  - Existing Rust kind: `rust-native` (academic algorithm; no specific C upstream being ported)
  - Existing non-C alternatives: —
  - Parallelism: `fastbloom` supports full concurrency via `portable-atomic`
  - SIMD: bit-manipulation auto-vectorizes; `fastbloom` documents this as its perf edge
  - Quadrant: ① (`fastbloom`) / ④ (others)
  - GPU-amenable: no — latency-bound point queries
  - Upstream license: public domain in spirit; `fastbloom` is `MIT OR Apache-2.0`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Adopt `fastbloom` for hot paths, fall back to `probabilistic-collections` for the feature-rich API. Used by ABySS, sourmash, many metagenomics tools.

- [x] **Cuckoo filter** — Bloom alternative with deletion + better locality.
  - Reference impl: `C++` · [efficient/cuckoofilter](https://github.com/efficient/cuckoofilter) · `Apache-2.0`
  - Existing Rust: [`cuckoofilter`](https://crates.io/crates/cuckoofilter) `0.5.0`; [`probabilistic-collections`](https://crates.io/crates/probabilistic-collections) `0.7.0`; [`autoscale_cuckoo_filter`](https://crates.io/crates/autoscale_cuckoo_filter) `0.5.21`
  - Existing Rust kind: `rust-native` (algorithm published by efficient/cuckoofilter but the Rust crates are independent impls, not code-ports)
  - Existing non-C alternatives: —
  - Parallelism: single-threaded today; concurrent variant is an open opportunity
  - SIMD: none explicit
  - Quadrant: ④
  - GPU-amenable: no — latency-bound point queries
  - Upstream license: `Apache-2.0`
  - Priority: `P1`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Strong fit for k-mer deduplication where deletion is needed (streaming metagenomics).

- [ ] **Compacted de Bruijn graph** — adjacency structure for assembly + pangenome.
  - Reference impl: `C++` · [GATB-bcalm](https://github.com/GATB/bcalm) · `MIT`
  - Existing Rust: [`debruijn`](https://crates.io/crates/debruijn) `0.3.4` ([`10XGenomics/rust-debruijn`](https://github.com/10XGenomics/rust-debruijn)); [`ggcat`](https://github.com/algbio/ggcat) (binary tool, not a published library on crates.io — install from source); [`rust-mdbg`](https://github.com/ekimb/rust-mdbg) (binary tool; companion library [`rust-seq2kminmers`](https://crates.io/crates/rust-seq2kminmers) `0.1.0`)
  - Existing Rust kind: `pure-port` (compaction in `ggcat` follows the bcalm/MEGAHIT algorithmic blueprint)
  - Existing non-C alternatives: —
  - Parallelism: explicit rayon-equivalent (`parallel-processor` in ggcat); 10x's `debruijn` is multi-threaded build
  - SIMD: ggcat uses `streaming-libdeflate-rs` and other SIMD-aware deps
  - Quadrant: ①
  - GPU-amenable: maybe — graph compaction is irregular, but k-mer counting prequel is SIMT-friendly
  - Upstream license: `MIT` (bcalm); ggcat is `MIT OR Apache-2.0`
  - Priority: `P1`
  - Layer: `A` (foundation — the graph type) + `B` (tool — `ggcat`-equivalent assembler)
  - Consumes primitives: `rsomics-kmer` (the ntHash crate), `rsomics-fm-index` indirectly
  - Notes: `ggcat` is production-grade and the right consolidation target. `debruijn` (10x) is older but still maintained. Note: the crates.io name `ggcat` is squatted by an unrelated crate (`ggcat = "0.0.1"`, clipboard tool); install algbio's ggcat from source. Cross-references [`02-genomics/assembly.md`](../02-genomics/assembly.md).
