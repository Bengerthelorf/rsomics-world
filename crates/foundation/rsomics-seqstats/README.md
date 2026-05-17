# rsomics-seqstats

Layer-A foundation primitives for sequence statistics, shared by the
`rsomics-*-stats` tools (`rsomics-fasta-stats`, `rsomics-fastq-stats`, …).
Library-only — no binary, no upstream to bench, perf-exempt by nature.

It is pure, format-agnostic math: feed it a `Vec<u64>` of record lengths and a
byte slice and it returns the seqkit-compatible length distribution and base
composition.

## API

- `LengthStats::new(Vec<u64>)` → `q1()` / `q2()` / `q3()` (length quartiles),
  `n50_l50()` → `(N50, L50)` where `L50` is seqkit's `N50_num`
  (unique-length-bucket count, not record count).
- `count_any_of(haystack, needles)` — deduped multi-byte count
  (GC bases, N runs, gap letters).
- `classify(sample)` / `SeqType` — seqkit's alphabet guess
  (DNA / RNA / Protein / Unlimit), `Serialize` with seqkit's exact JSON
  names. seqkit guesses from the **first record only**, scanning at most
  `DEFAULT_ALPHABET_GUESS_LEN` bytes — callers pass the first record's
  sequence prefix, never a cross-record accumulation.

## Origin

This crate is a Rust port of the length-distribution math in
[`shenwei356/bio`](https://github.com/shenwei356/bio) `util/length-stats.go`
and the alphabet guess in [`seqkit`](https://github.com/shenwei356/seqkit)
`seq.GuessAlphabet`, so `rsomics-*-stats --all --tabular` byte-agrees with
`seqkit stats -a -T`. seqkit computes L50 over unique-length buckets rather
than records; that quirk is reproduced deliberately for compatibility.

Both upstreams are MIT, so the source was read and is cited directly (the
clean-room rule applies only to GPL upstreams).

License: MIT OR Apache-2.0.
Upstream credit: shenwei356/bio, shenwei356/seqkit (MIT).
