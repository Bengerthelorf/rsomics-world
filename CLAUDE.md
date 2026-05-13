# rsomics-world — Autopilot Operating Manual

You are running in **autopilot mode** on the `rsomics-world` repository at `/Volumes/Zane's HDD/Documents/rsomics-world/`. This file is loaded every session — it is the constitution. The current Phase's prompt is the law of the current period.

## Mission

Plan and build a **monorepo cargo workspace** of `rsomics-<name>` single-binary CLI tools — "many bcmrs" — that displace the C/Python/R-era bioinformatics toolchain. The motivation is concrete: most upstream tools are single-threaded, memory-inefficient, written in 2005-era C or pure R, and waste modern multicore + SIMD + GPU resources.

Quality bar per tool: **fail-fast, fail-loud, single binary, easy to install (`cargo install rsomics-<name>`), zero defensive programming, self-explaining code**.

## Architecture (2 layers, monorepo)

```
rsomics-world/
├── Cargo.toml              ← workspace manifest
├── CLAUDE.md               ← this file
├── CONVENTIONS.md          ← repo rules (autopilot may evolve, with care)
├── ROADMAP.md / TODO.md    ← planning catalog (kept current)
├── README.md               ← user-facing project intro
├── docs/                   ← planning per module (00-overview … 09-workflow-utility)
├── crates/
│   ├── foundation/         ← Layer A: library-only, depended on by 2+ tools
│   │   ├── rsomics-common/     (errors, CLI scaffold, runner, json/progress/exit)
│   │   ├── rsomics-intervals/  (BED algebra)
│   │   ├── rsomics-kmer/       (hashing, counting, sketches)
│   │   ├── rsomics-fm-index/   (BWT, suffix array, locate)
│   │   ├── rsomics-align-core/ (SW/NW/WFA/banded)
│   │   └── rsomics-stats/      (GLM, FDR, p-values)
│   └── tools/              ← Layer B: each is one installable binary
│       ├── rsomics-fastp/
│       ├── rsomics-bam/    (subcommands: view, sort, index, markdup, ...)
│       └── ...
└── .autopilot/             ← persistence (gitignored)
    ├── sessions/           ← per-session log
    ├── gates/              ← gate status reports awaiting review
    ├── needs-review/       ← external claims that failed verification gates
    └── state/              ← long-running cache (classified entries, etc.)
```

### Dependency rules (one-way, enforced)

- B → A → external. **A never depends on B. B never depends on B directly** — share via A.
- A is **library-only** (no `[[bin]]`). A crate with a binary is by definition B.
- Promote internal module to A only when **2+ B crates** need it (YAGNI; do not pre-abstract).

### External dependency 4-quadrant classification

For every adopted external crate, classify and document in its TODO entry:

| Quadrant | What it is | Acceptable for | Note |
|---|---|---|---|
| ① Pure Rust + explicit parallelism (rayon / std::simd) | E.g. noodles, rayon, polars, needletail | Adopt freely | This is the goal state |
| ② FFI wrapper over C (`cc`, `bindgen`, `*-sys`) | E.g. rust-htslib, minimap2-rs | Adopt for speed, but **document clearly: "FFI wrapper, not Rust rewrite"** | C upstream stays the perf reference; future P2 may pure-Rust replace |
| ③ Pure Rust but single-threaded in hot path | (varies — must verify per case) | **Avoid** in hot paths | Inherits the disease we're curing |
| ④ Pure Rust, non-hot, edge utility | clap, serde, anyhow, regex | Adopt freely | Doesn't matter |

When adopting, the TODO entry's `Existing Rust kind` field records ①②③④ explicitly.

## Operating mode

- **No clarifying questions.** Decide, commit, document the reason in the commit body.
- **Don't stop unless a real decision is needed.** Autopilot means autopilot. Inside a phase, when the directive is clear (e.g. "implement rsomics-fastp's preprocessing pipeline"), keep going through the next logical sub-task without surfacing for confirmation. **Halt only when**: (a) you hit one of the gate-defined stop conditions or halt triggers, (b) a genuine decision arises that materially changes scope or requires user judgment (not a stylistic preference, not a "good stopping point"), or (c) you finish the phase. Session-end courtesy halts ("opening landed, ready for next session") are not allowed — if there is more work scoped within the current phase / directive, do it. If the user wanted you to stop, they would say so.
- **Solo, no PRs.** Direct commits to `main`. Conventional prefixes:
  - `docs(<module>):` for catalog content (under `docs/`)
  - `feat(<crate>):` for new crate or feature
  - `fix(<crate>):` for bug fix
  - `refactor(<crate>):` for restructuring
  - `test(<crate>):` for tests
  - `chore:` for housekeeping
  - `ci:` for `.github/workflows/`
- **No `Co-Authored-By` trailer. Ever.**
- **One concern per commit. Fine-grained, push every 3-5 commits.** After every push: `gh run list --branch main --limit 3` and wait for green before tagging.
- **Tags are user decisions** — never tag, never `gh release create`. Surface to user when a phase deserves a release.

## Gate-based checkpoints

Autopilot is autonomous **within a phase**. Between phases, halt and let the user review.

Each phase ends by writing `.autopilot/gates/gate-<N>-<date>.md` summarizing:
1. What changed (file list + commit shas).
2. **Decisions taken without asking (with reasoning).** List every autonomous choice — version pins, runner-label substitutions, default values picked, schema additions, scope expansions, naming choices — even when the choice feels like a neutral substitution or an improvement. Distinguish direct user policy ("user said untrack ROADMAP.md") from autopilot autonomy ("substituted retired macos-13 with macos-15-intel"). User policy doesn't go in this section; autonomy does.
3. Needs-review backlog (count + categories).
4. What the next phase requires the user to bless.
5. Recommended next phase prompt.

Then **stop**. Next session: read the most recent un-approved gate file; if the user has appended `approved: yes`, proceed. Otherwise wait or revise per their comments.

## Verification gates — non-negotiable, applied to EVERY external claim

For each crate / repo / version mentioned in a TODO entry, run all five before commit:

1. **Existence**: `cargo search <name> --limit 3`. Empty → claim stale → log + skip.
2. **Aliveness**: `gh repo view <owner>/<repo> --json updatedAt,description,isArchived,defaultBranchRef`. Archived or > 18 months stale → downgrade or skip.
3. **Binary shape**: `WebFetch https://raw.githubusercontent.com/<owner>/<repo>/<branch>/Cargo.toml`. Look for `[[bin]]`, `default-run`. Record `ships-binary` or `library-only`.
4. **Substance**: `WebFetch` `src/main.rs` or `src/lib.rs`. < ~500 LOC or `unimplemented!()` → stub → downgrade.
5. **Perf-class** (the new gate): in Cargo.toml grep for `cc`, `bindgen`, `*-sys`; check `build.rs` exists. In lib/main.rs grep for `rayon`, `par_iter`, `std::simd`, `packed_simd`, `target_feature`. Fill quadrant ①②③④ in entry.

**Any fetch fails** (404, network, blocked) → append to `.autopilot/needs-review/external-<date>.md` and skip this entry. **NEVER invent a crate name, URL, or version.**

## Code style (enforced in EVERY commit)

- **Fail-fast, fail-loud.** `Result` propagation, no `unwrap_or_default()` to hide errors, no defensive `if let Ok()` swallows. Errors travel to top-level main; main exits non-zero with stderr message. **Wrong output is worse than crash in bioinformatics — bail rather than ship a wrong VCF.**
- **No defensive programming.** Don't validate invariants the type system enforces; don't double-check things already checked at the boundary; don't add retries / fallbacks for cases that can't happen.
- **Comments are persistent.** No `// Phase 1 / Phase 2 / Phase A / Phase B`, no `// TODO: do later`, no `// removed X`, no `// based on the audit`. If a comment references a temporal state, it's wrong. Comments describe **invariants and non-obvious WHY**, not timelines.
- **No what-comments.** If a comment restates code, delete it. Names carry meaning.
- **Comment heuristic**: ask "would this comment still make sense if a fresh reader saw the code two years from now with no project history?" If no, delete.
- **Diagrams**: D2 or Mermaid only. Never ASCII art. (Pre-existing ASCII in `ROADMAP.md` may stay until that file is edited; if edited, convert.)
- **No meta-docs added.** CHANGELOG.md / RELEASING.md / ARCHITECTURE.md do NOT belong in this repo. The release notes live on the GitHub Release. Architecture lives in `CONVENTIONS.md` + `docs/00-overview/`.
- **Unwrap discipline**: prod code uses `unwrap()` only when the invariant is statically obvious (e.g., `slice.first().unwrap()` after a non-empty check at the type level). Tests can unwrap freely.

## Compatibility + perf — the rsomics-* contract

Every Layer B crate that ports a C/Python/R tool MUST have:

- `tests/golden/` — small (KB-MB) reference inputs.
- `tests/compat.rs` — runs the rsomics-* binary, runs the upstream binary, byte-or-field-level diff. **Fail loud on mismatch.**
- `benches/<tool>.rs` — criterion bench against the upstream binary on a non-trivial input. CI tracks regression > 10% on hot path → fail.
- README "Origin" section (license clean-room methodology — see below).

No crate ships without these three. A crate with only `tests/` and no `compat.rs` is an unfinished implementation.

## License + clean-room methodology

Upstream tools have varied licenses. Default for our crates: **MIT OR Apache-2.0**.

For **GPL upstream** (bowtie2, HISAT2, original BWA, MEGAHIT, MMseqs2, Trimmomatic, Subread/featureCounts, SPAdes, IQ-TREE, etc.):

- Our Rust port can still be MIT/Apache-2.0 **only if** the implementation is clean-room: based on the published paper + format spec + black-box behavior observation, NOT on reading the GPL source.
- Each such crate's README MUST contain:

  ```markdown
  ## Origin

  This crate is an independent Rust reimplementation of `<upstream>` based on:
  - The published method (cite paper + DOI)
  - The public file-format spec
  - Black-box behavior testing against the upstream binary

  No source code from the GPL upstream was used as reference during
  implementation. Test fixtures are independently generated or sourced
  from public benchmark datasets (HG002 / 1000 Genomes / etc.).

  License: MIT OR Apache-2.0.
  Upstream credit: <upstream> <link> (<their license>).
  ```

- If clean-room cannot be confirmed for an entry, the entry is downgraded to: defer or write as GPL crate. Flag and surface — do not assume.

For **FFI-wrapper adopt** (rust-htslib, minimap2-rs etc.): record in entry as **Quadrant ②**, credit upstream prominently, document license inheritance.

## Cross-platform target

- **First-class**: `x86_64-unknown-linux-gnu`, `aarch64-unknown-linux-gnu`, `x86_64-apple-darwin`, `aarch64-apple-darwin`. CI runs all four. Bench runs Linux x86_64 + aarch64.
- **Not supported**: Windows (skip in CI), other OSes. If a tool's algorithm hinges on Linux-only syscalls (e.g., `splice`), gate behind `#[cfg(target_os = "linux")]` with a portable fallback.

## CI policy (GitHub Actions)

Workflows in `.github/workflows/`:

- `ci.yml`: fmt + clippy + test matrix (4 first-class targets × stable + MSRV).
- `bench.yml`: criterion benches on Linux x86_64 + aarch64, manual trigger. Regression detection.
- `release.yml`: triggered by tag push; cross-compiles binaries, uploads to GitHub Release.

**CI is the truth.** Local Mac smoke is supplementary. After every push, `gh run list --branch main --limit 3` and wait for green before any tag.

## Disk + environment hygiene (mini_m2 + autopilot)

- **mini_m2's boot disk is small.** All cargo target, all build artifacts, all caches → HDD at `/Volumes/Zane's HDD/`.
- **Session bootstrap (every shell)**:
  ```bash
  export CARGO_TARGET_DIR="/Volumes/Zane's HDD/cargo-target/rsomics-world"
  export CARGO_HOME="/Volumes/Zane's HDD/cargo-home"
  cd "/Volumes/Zane's HDD/Documents/rsomics-world"
  ```
- Run this before every `cargo` invocation. If `df -h /` shows boot disk filling, **halt** and surface.
- `.gitignore`: `/target`, `.autopilot/`, `Cargo.lock` only for binary crates (workspace root keeps the lockfile).

## Test data tiers (CPU/disk-aware)

| Tier | Where | Size | Use |
|---|---|---|---|
| 1 | inside git (`tests/golden/`) | < 100 MB total | Small synthetic FASTQ/BAM/VCF for unit tests |
| 2 | GitHub Actions runner (downloaded) | < 5 GB ephemeral | chr20/22-subset HG002, public 1000G subsets — integration tests |
| 3 | mini_m2 HDD (`/Volumes/Zane's HDD/rsomics-fixtures/`) | 50-200 GB | One real WGS sample, a few RNA-seq, smoke benchmarks |
| 4 | `ssh 4090:/data3/rsomics-fixtures/` | up to ~1 TB | Multi-sample cohort benches, large-scale perf validation |

Tier 1 in git. Tier 2 downloaded on demand (manifest in repo lists URLs + sha256). Tiers 3-4 manually provisioned; autopilot references via `BCMR_BENCH_DATA` env var.

## GPU strategy

**CPU first, always.** Each TODO entry's Notes records `GPU-amenable: yes/no/maybe` with one-sentence reason, but GPU implementation is gated on:
1. CPU implementation complete + passing compat tests.
2. CPU perf baseline measured.
3. GPU is a measurable win (not a religious win).

GPU testing uses `ssh 4090` (`/data3` for fixtures, `nvidia-smi` available). Framework: `candle` for inference, `burn` if cross-backend is needed.

## Autopilot memory hygiene

- Write `project_*` memory entries for non-obvious discoveries (e.g., "crate X claims pure-Rust but `build.rs` compiles 2000 lines of C — Quadrant ②").
- Don't write transient state (current task, what file I'm editing) to memory — that's session log territory.
- At each gate transition, scan `~/.claude/projects/.../memory/MEMORY.md` for stale entries and trim. Stale = contradicted by current code, or referencing a phase that's over.

## Toolbox (quick reference)

```bash
# verification
cargo search <name> --limit 3
gh repo view <owner>/<repo> --json updatedAt,description,isArchived
WebFetch https://raw.githubusercontent.com/<owner>/<repo>/<branch>/Cargo.toml
WebFetch https://raw.githubusercontent.com/<owner>/<repo>/<branch>/src/{main,lib}.rs

# cargo on mini_m2 — always with HDD target
export CARGO_TARGET_DIR="/Volumes/Zane's HDD/cargo-target/rsomics-world"
cargo build --release
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all

# git
git log --oneline -10
git status
gh run list --branch main --limit 3

# remote bench
ssh 4090 "cd /data3/rsomics-fixtures && ./run-bench.sh"
```

## FreshEye discipline (phase-conditional)

Autopilot cannot audit its own systematic errors. Multi-model / multi-session review catches what self-audit misses. Apply at intensity matching each phase's blast radius.

### Levels

**L0 — no FreshEye** — purely mechanical work (formatting, commit hygiene). Cost > value.

**L1 — sampling spot-check** — every N outputs, dispatch a fresh subagent (different model where possible) to re-audit a random ~10% sample. Mismatches → `.autopilot/needs-review/<topic>-<date>.md` for batched user review. Default for **Phase 1 catalog work** and **non-load-bearing doc edits**.

**L2 — full per-output review** — every significant output (a new public API in a foundation crate, a Layer A crate's `lib.rs`, a non-trivial commit) reviewed end-to-end by a fresh subagent before commit. Used in **Phase 2-3** for foundation crates whose API is downstream-critical.

**L3 — multi-axis parallel review** — dispatch 2-3 fresh subagents on different lenses (reuse / quality / efficiency / correctness) before commit. Used in **Phase 4 gate (killer binary)** and on any Layer A crate that ≥ 3 tools will depend on.

### How to dispatch

**Anthropic models (abundant, primary)** — use the Agent tool. Sonnet quota is unconstrained for this project (explicit user policy — prefer Sonnet for routine FreshEye work):

```
Agent(
  subagent_type: "general-purpose",
  model: "sonnet" | "opus" | "haiku",
  description: "FreshEye <level>: <topic>",
  prompt: "<self-contained brief: what to audit, what to check for, output format>"
)
```

**Codex (OpenAI GPT-5.x, quota-constrained)** — dispatch as a subagent via the Agent tool, but put `--effort medium` on the first line of the prompt. **Never `high` / `xhigh`** — quota expensive:

```
Agent(
  subagent_type: "codex:codex-rescue",
  description: "FreshEye <level> (codex): <topic>",
  prompt: "--effort medium\n\n<self-contained brief>"
)
```

First use of the session must precede with `Skill(codex:setup)` to confirm CLI ready. Spend Codex quota on: Phase 4 killer-binary L3 reviews, L2-triggered foundation-crate API reviews, cross-family second opinions. **Phase 1 catalog work does not use Codex** — Sonnet ↔ Opus internal review is sufficient.

**Gemini (Google Gemini 3.x, quota abundant)** — Gemini is a **Skill**, not a subagent. Dispatch via the Skill tool. The first line of `args` selects the model — `gemini-3.1-pro` for heavy review, flash variants for light review:

```
Skill(
  skill: "gemini:rescue",
  args: "Use gemini-3.1-pro. <self-contained brief>"
)
```

First use of the session must precede with `Skill(gemini:setup)` to confirm CLI ready. **Gemini quota is fine to spend freely**, particularly well-suited for tests / fuzz / Cargo.toml audits — axes that Anthropic models tend to under-emphasise.

**Cross-family selection rules**:

- Same-output **never** self-audited — Sonnet ↔ Opus internal pair-review is the default.
- Adding a non-Anthropic axis: `tests/` / non-production code / Cargo.toml → Gemini; `src/` production code fresh-eye → Codex medium.
- All three families: Phase 4 killer-binary L3 — Sonnet + Opus + Codex medium + Gemini-pro for four-model differentiated perspectives.

**Setup failure fallback**: if `codex:setup` / `gemini:setup` is unavailable in this session, FreshEye degrades to Anthropic-only (Sonnet ↔ Opus) — phase work continues but the gate report records `cross-family review unavailable this session, fallback to Anthropic-only` so review-strength downgrades surface to the user.

**Quota tracking**: autopilot maintains a running tally in `.autopilot/state/freshseye-budget.toml` of `codex-runs` and `gemini-runs` for the session, so the user can judge whether Phase 4's L3 budget still has room.

Never run FreshEye against output produced in the same context — `fresh` means fresh context.

### Mandatory triggers (any phase, regardless of declared level)

- About to commit a `pub` item in a foundation crate → L2 minimum.
- About to publish to crates.io → L3.
- Gate report drafted from > 200 individual decisions → L1-sample those decisions retroactively.
- Real-world test contradicts unit tests (especially platform-quirky details) → dispatch debug-focused FreshEye agent.
- The instruction `external_advice_must_question_assumption` applies recursively: a FreshEye finding is itself external advice and must be sanity-checked against deployment reality before action.

### Real-world testing — FreshEye's twin

Local CI + unit tests are necessary but not sufficient. For anything touching cross-platform / SSH / GPU / FS quirks, exercise on:

- mini_m2 (`/Volumes/Zane's HDD/`) — macOS + local FS
- `ssh 4090` (`/data3` fixtures, GPU) — Linux + perf + GPU
- GitHub Actions matrix — broad targets, no GPU

If a test passes on one and fails on another, the test is fragile, not the code.

## Universal stop conditions (any phase)

- About to commit a claim you couldn't verify → STOP, log to needs-review.
- About to invent a crate name / URL / version → STOP, log to needs-review.
- Three external signals diverge for one tool (crates.io / GitHub / source disagree) → STOP, document.
- More than 5 unrelated files in one commit → wrong scoping, restart.
- More than 3 iterations producing no meaningful diff → easy wins exhausted, gate out.
- Disk usage on `/` > 80% → STOP, audit cargo target placement.
- `needs-review/` accumulates > 20 entries → too many unverifiables, halt for batch review.
- Architecture (2-layer / 4-quadrant / monorepo) looks wrong for the task at hand → never restructure unilaterally; surface to user.
