# rsomics-fastq-correct

K-mer-spectrum substitution-error correction for Illumina FASTQ reads
(SE and PE) — an independent Rust port of BFC.

```bash
cargo install rsomics-fastq-correct
```

## Scope

This crate is the **read-error-correction** partition of the FASTQ
preprocessing surface (one operation = one crate). It counts every input
k-mer, derives a trusted-k-mer set from the coverage spectrum, and
rewrites untrusted bases along the highest-support correction path. It is
a standalone pre-assembly / pre-variant-calling pipeline stage.

| Operation | Crate |
|---|---|
| K-mer-spectrum error correction | **rsomics-fastq-correct** ← here |
| Adapter / poly-X / fixed trim | rsomics-fastq-trim |
| Quality + length filter | rsomics-fastq-filter |
| Inline-UMI extract + stamp | rsomics-fastq-umi |
| Exact / near dedup | rsomics-fastq-dedup |

## Behaviour

Faithful to BFC 0.1 (`bfc_opt_init` defaults `-k 33`, `-c 3`, `-w 10`,
`-q 20`; the fixed penalty weights `w_ec=1 w_ec_high=7 w_absent=3
w_absent_high=1`, `max_path_diff=15`, `max_heap=100`, `max_end_ext=5`):

- One count pass over the whole input builds the trusted-k-mer table.
- Per read: the longest run of solid k-mers (count ≥ `--min-cov`) anchors
  a best-first correction walk (a min-heap on the BFC weighted penalty),
  run forward then on the reverse complement; the two directional results
  are merged by BFC's agreement rule.
- A read with no solid k-mer is not given up on immediately: BFC's
  greedy single-substitution probe (`bfc_ec_greedy_k`) scans successive
  first-k-mers for one confidently-supported base fix and re-anchors on
  success.
- Reads with > 5 % ambiguous bases, an uncorrectable `N`, too many
  dead-ends, or no solid k-mer even after the greedy probe are
  uncorrectable: passed through unchanged, or dropped with
  `--discard-uncorrectable`.
- Corrected bases are lower-cased; with qualities, a corrected base's
  quality byte encodes the original base (BFC's `34 + ob` convention).

## Usage

```bash
rsomics-fastq-correct -i reads.fq.gz -o corrected.fq.gz          # SE, k=33
rsomics-fastq-correct -i r1.fq.gz -I r2.fq.gz \
                       -o c1.fq.gz -O c2.fq.gz -c 5 -D --json     # PE, strict
```

## Origin

This crate is an independent Rust reimplementation informed by reading
the permissively-licensed upstream source (allowed and cited):

- **BFC** — `bfc.c`, `correct.c`, `kmer.h`, `bfc.h` read at the `master`
  tag ([lh3/bfc](https://github.com/lh3/bfc), MIT). The trusted-k-mer
  decision, the `bfc_ec1dir` best-first penalty search, the bidirectional
  `bfc_ec1` orchestration, the `bfc_kmer_hash` canonical key, and every
  numeric default/penalty weight are ported from this source. Paper:
  Li 2015, *Bioinformatics*, doi:10.1093/bioinformatics/btv592.
- **Lighter** ([mourisl/Lighter](https://github.com/mourisl/Lighter),
  GPL-3.0) is a secondary *behavioural* reference only — no GPL source
  was read; clean-room (paper + black-box) applies to it.

### Implementation decisions (documented divergences)

1. **Count table.** BFC's counter is its own open-addressing
   `bfc_ch` hash table sized by `-b`/`-H`/`-s`. We key a plain map by
   BFC's exact `bfc_kmer_hash` value, so the *trusted-k-mer counts*
   match BFC's, but the table's collision behaviour differs. BFC's
   bloom-sizing / hash dump-restore / chunk / verbosity flags
   (`-b -H -s -d -r -L -V`) are specific to `bfc_ch` and are
   intentionally not surfaced — they have no semantics in this
   structure. (`-s` genome-size auto-thresholding is a distinct
   sub-mode; the explicit `--min-cov` path is the primary surface.)
2. **Compatibility is semantic, not byte-identical.** BFC's correction
   tie-breaks depend on its specific counter's collision profile, so a
   byte-exact diff against `bfc` is not well-defined across counter
   implementations. `tests/compat.rs` asserts *correction-outcome
   equivalence* on golden reads whose k-mer coverage makes the corrected
   base implementation-independent, version-gated against `bfc`. This is
   a documented justified divergence (cf. `rsomics-fastq-split`'s
   exact-vs-heuristic `--split N`). The perfgate is still strict
   `> 1.0×` throughput vs `bfc`.

License: MIT OR Apache-2.0.
Upstream credit: [BFC](https://github.com/lh3/bfc) (MIT),
[Lighter](https://github.com/mourisl/Lighter) (GPL-3.0).
