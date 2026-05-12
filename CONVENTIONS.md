# Conventions

Rules for contributing to this planning repo and for the Rust crates it will
eventually spawn.

## Doc format

Every module sub-doc (the topic files under `docs/0X-*/`) follows this skeleton:

```markdown
# <Topic title>

> One-sentence elevator description.

## Scope

What problems this topic covers; what it does *not* cover (point to the
neighboring topic).

## Design notes

Algorithmic considerations, where Rust helps, where it doesn't. Brief — this is
not a research paper. 3–6 bullets is typical.

## TODO

A flat checklist. Each entry has the structure below.
```

### TODO entry format

```markdown
- [ ] **<Canonical tool name>** — <one-line purpose>
  - Reference impl: <language> · <repo URL or "—"> · <license>
  - Existing Rust: <crate name + URL, or "none">
  - Existing non-C alternatives: <Zig / Go / C++ rewrites worth knowing, or "—">
  - Priority: P0 / P1 / P2  (P0 = must have, P1 = high value, P2 = nice to have)
  - Notes: <SIMD-critical? GPU? Already production-ready in another language?
    Should we wrap vs rewrite?>
```

**Rule:** every tool gets an entry, *even if a mature Rust implementation
already exists*. Mark such entries with `[x]` and explain in the notes whether
we plan to (a) adopt as-is, (b) extend, or (c) leave alone.

## Crate naming (future code)

When we start writing crates:

- Top-level umbrella crate: `rsomics`
- Domain crates: `rsomics-<domain>` (e.g. `rsomics-io`, `rsomics-align`,
  `rsomics-vcf`, `rsomics-scrna`).
- Binary tools shipped alongside libraries get the same prefix:
  `rsomics-bwa`, `rsomics-samtools`.
- Each crate lives in its **own repository** under the same GitHub org/user,
  *not* as a folder in this repo. This planning repo stays docs-only.

## License

- Docs in this repo: **CC BY 4.0**.
- Each future crate: dual **MIT OR Apache-2.0** (standard Rust ecosystem
  pattern), unless an upstream tool we are deriving from forces otherwise (GPL
  derivatives must be flagged in the relevant module doc).

## Linking between docs

Use relative paths: `[topic](../02-genomics/alignment-short-read.md)`.
Cross-module references are encouraged — many tools straddle categories.

## Adding a new topic

1. Pick the module it belongs in. If it does not fit, raise an issue rather
   than inventing a 10th top-level module.
2. Add the file using the doc-format skeleton above.
3. Update the module's `README.md` index.
4. Add the entries to `TODO.md` (or run the aggregator once we write one).
