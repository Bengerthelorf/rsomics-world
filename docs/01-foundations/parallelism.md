# Parallelism

> Threading patterns, async I/O, and GPU offload — how rsomics crates exploit
> modern hardware.

## Scope

Patterns and primitives rather than tools. `rayon` for data-parallel work,
`tokio` (rarely) for I/O-bound concurrency, `crossbeam` channels for
pipeline stages, and Rust GPU stacks (`candle`, `burn`, `wgpu`, `cust`).
Excludes specific GPU-accelerated tools (those live in their domain module —
e.g. GPU AlphaFold inference under [`07-proteomics-structure/`](../07-proteomics-structure/)).
The benchmarking framework lives in
[`00-overview/benchmarking-plan.md`](../00-overview/benchmarking-plan.md).

## Design notes

- Most bioinformatics workloads are embarrassingly parallel over records.
  `rayon::par_iter()` reaches near-linear scaling on 64-core boxes when
  records are independent (FASTQ → trimmed FASTQ, BAM → tagged BAM, VCF
  → annotated VCF). This is the default tool.
- Where Rust beats htslib-era C: **scaling**. Most htslib commands
  parallelise only BGZF decompression (a single `-@`/`--threads` knob),
  not the work that follows. As of 2024, htslib still only multi-threads
  BAM compression on the *write* side, not BAM decompression on read.
  Rust pipelines can fan out the whole record processing stage.
- Where Rust struggles: **I/O-bound steady-state** (e.g. an aligner
  feeding from a remote BAM). `tokio` adds complexity that rarely pays
  off vs. a tuned thread pool with blocking reads. Use async only when
  cloud streaming dominates wall-time.
- GPU offload is a second-stage concern (per
  [`ROADMAP.md`](../../ROADMAP.md): "GPU is Phase 4+"). The Rust GPU
  ecosystem (candle, burn-wgpu, cust) is improving fast — adoptable for
  deep-learning-based callers (DeepVariant-class) and AlphaFold inference,
  but not for classical aligners where SIMD CPU is still cheaper per dollar.
- Document the threading model in every crate's README. A user should
  know whether to set `RAYON_NUM_THREADS`, `--threads`, both, or neither.

## How big tools thread today

| Tool       | Threading model | Notes |
|------------|-----------------|-------|
| `samtools` | pthreads, BGZF + per-command pool | Scales to ~16 threads on most commands, sub-linear past that |
| `bcftools` | pthreads | Mostly single-threaded outside compression |
| `bwa-mem`  | pthreads, per-read | Near-linear to ~32 threads, then memory-bandwidth-bound |
| `bwa-mem2` | pthreads + AVX-512 SIMD | Beats bwa-mem 2-3× single-thread; threading model same |
| `minimap2` | pthreads + SSE/AVX | Near-linear to ~32 threads |
| `GATK HC`  | Java Spark / native threads | Spark mode adds JVM/network overhead |
| `STAR`     | pthreads + huge shared index | Memory-bound past 8 threads |

## TODO

- [x] **`rayon`** — work-stealing data parallelism.
  - Reference impl: — (Rust-native; analogous to Intel TBB / OpenMP)
  - Existing Rust: [`rayon`](https://github.com/rayon-rs/rayon) (production)
  - Existing non-C alternatives: TBB (C++), OpenMP, Go goroutines
  - Priority: `P0`
  - Notes: Adopt as the default. Convention: `par_iter()` is the first
    tool reached for; document `RAYON_NUM_THREADS` interaction with any
    explicit `--threads` flag.

- [x] **`crossbeam-channel`** — bounded MPMC channels for pipeline stages.
  - Reference impl: — (Rust-native; analogous to Go channels)
  - Existing Rust: [`crossbeam-channel`](https://crates.io/crates/crossbeam-channel)
  - Existing non-C alternatives: Go channels; C++ `tbb::concurrent_bounded_queue`
  - Priority: `P0`
  - Notes: Adopt for reader → worker → writer pipelines (the standard
    bioinformatics streaming pattern). Avoid `std::sync::mpsc` (slower,
    SPMC only).

- [~] **`tokio` async I/O** — for cloud/network-bound workloads.
  - Reference impl: — (Rust-native; analogous to Go runtime)
  - Existing Rust: [`tokio`](https://github.com/tokio-rs/tokio) (production)
  - Existing non-C alternatives: Boost.Asio (C++), Go runtime
  - Priority: `P1`
  - Notes: noodles ships async variants of most readers behind a feature
    flag. Use sparingly — only when the workload is genuinely
    I/O-bound (htsget streaming, S3 / GCS BAM access). CPU-bound code
    does *not* benefit from async.

- [~] **`candle`** — pure-Rust deep-learning framework.
  - Reference impl: `Python/C++` · [pytorch/pytorch](https://github.com/pytorch/pytorch) · `BSD-3-Clause`
  - Existing Rust: [`candle`](https://github.com/huggingface/candle) (CPU, CUDA, Metal)
  - Existing non-C alternatives: PyTorch, JAX, MLX (Swift)
  - Priority: `P1`
  - Notes: First-choice for DL inference (DeepVariant-class callers,
    AlphaFold). Smaller and faster-loading than PyTorch but covers fewer
    ops. Track `candle-transformers` for ESM/AlphaFold-relevant layers.

- [~] **`burn`** — alternative pure-Rust DL framework with backend
  swapping.
  - Reference impl: `Python/C++` · [pytorch/pytorch](https://github.com/pytorch/pytorch) · `BSD-3-Clause`
  - Existing Rust: [`burn`](https://github.com/tracel-ai/burn) (`burn-wgpu`, `burn-candle`, `burn-ndarray`, `burn-tch` backends)
  - Existing non-C alternatives: PyTorch, TensorFlow
  - Priority: `P1`
  - Notes: More framework-y than candle (full training loop, autodiff,
    dashboard). The wgpu backend means cross-vendor GPU (NVIDIA / AMD /
    Apple / WebGPU) without CUDA dependency — useful for portable
    inference binaries.

- [~] **`wgpu`** — portable GPU compute via WebGPU shading language.
  - Reference impl: — (browser standard, native via Vulkan/Metal/D3D12)
  - Existing Rust: [`wgpu`](https://github.com/gfx-rs/wgpu) (production)
  - Existing non-C alternatives: CUDA (NVIDIA only), Vulkan compute
  - Priority: `P2`
  - Notes: Low-level GPU access for custom kernels (Smith-Waterman,
    rolling hash) when DL frameworks are overkill. Adopt indirectly via
    `burn-wgpu`; direct use is a Phase 4 concern.

- [ ] **CUDA via `cust`** — direct NVIDIA GPU compute.
  - Reference impl: `C++/CUDA` · [NVIDIA/cccl](https://github.com/NVIDIA/cccl) · `Apache-2.0`
  - Existing Rust: [`cust`](https://github.com/Rust-GPU/Rust-CUDA);
    [`cudarc`](https://github.com/coreylowman/cudarc)
  - Existing non-C alternatives: CUDA C++, JAX
  - Priority: `P2`
  - Notes: Only worth it when a tool ships a CUDA inner kernel (e.g.
    GPU BWA, GPU minimap2 forks). For DL inference, candle/burn handle
    CUDA without raw `cust`.

- [ ] **SIMD portable abstractions** — `std::simd` and the `core_arch`
  intrinsics.
  - Reference impl: `C++` · intrinsics + `std::experimental::simd`
  - Existing Rust: [`std::simd`](https://doc.rust-lang.org/std/simd/index.html) (nightly);
    [`wide`](https://crates.io/crates/wide);
    [`pulp`](https://crates.io/crates/pulp) (runtime dispatch)
  - Existing non-C alternatives: ISPC; Highway (C++)
  - Priority: `P0`
  - Notes: Per
    [`00-overview/design-principles.md`](../00-overview/design-principles.md)
    rule 3, every hot kernel has a scalar fallback and a SIMD path. `pulp`
    handles runtime feature detection nicely. Avoid AVX-512-only code
    paths in shipped binaries (still poorly supported on consumer CPUs).
