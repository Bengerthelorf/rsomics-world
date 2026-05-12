# <Topic title>

> One-sentence elevator description of what this topic covers.

## Scope

What this topic includes; what it explicitly excludes (and where that work
lives instead).

## Design notes

3–6 bullets. Algorithmic considerations, where Rust shines, where it
struggles, and which existing crates we expect to lean on or extend.

- Bullet one.
- Bullet two.
- Bullet three.

## TODO

Each entry uses the format below. Every canonical tool gets an entry, even
if a mature Rust port already exists — in that case mark `[x]` and explain
in Notes whether we adopt, extend, or leave it alone.

- [ ] **`<canonical-tool-name>`** — <one-line purpose>.
  - Reference impl: `<language>` · [`<repo>`](<URL>) · `<license>`
  - Existing Rust: `<crate or "none">` (link, status)
  - Existing non-C alternatives: `<Zig / Go / Julia / modern C++ rewrites>` or `—`
  - Priority: `P0` / `P1` / `P2`
  - Notes: <SIMD? GPU? FFI-only Rust binding so far? Should we wrap, rewrite,
    or skip?>

Example filled entry:

- [~] **`bwa-mem`** — Burrows-Wheeler short-read aligner.
  - Reference impl: `C` · [lh3/bwa](https://github.com/lh3/bwa) · `MIT`
  - Existing Rust: none pure-Rust; partial: [`bwa-mem2-rs`](https://github.com/) FFI only
  - Existing non-C alternatives: `bwa-mem2` (C++/SIMD rewrite, same author)
  - Priority: `P0`
  - Notes: Inner SW kernel benefits hugely from SIMD. Start with FFI wrapper,
    plan pure-Rust port after Phase 1 IO stabilizes. Compare against
    `bwa-mem2` not `bwa` 0.7.17 for fairness.
