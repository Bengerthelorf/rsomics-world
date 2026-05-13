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
- GPU offload is a second-stage concern: ramp begins after Phase 3. The
  Rust GPU ecosystem (candle, burn-wgpu, cust) is improving fast —
  adoptable for deep-learning-based callers (DeepVariant-class) and
  AlphaFold inference, but not for classical aligners where SIMD CPU is
  still cheaper per dollar.
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
  - Existing Rust: [`rayon`](https://github.com/rayon-rs/rayon) `1.12.0`
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: TBB (C++), OpenMP, Go goroutines
  - Parallelism: this **is** the parallelism crate (work-stealing thread pool, parallel iterators)
  - SIMD: none (this is the scheduler; SIMD lives in the kernels)
  - GPU-amenable: n/a — CPU scheduler
  - Upstream license: `MIT OR Apache-2.0`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Quadrant ①. Adopt as the default. Convention: `par_iter()` is the first tool reached for; document `RAYON_NUM_THREADS` interaction with any explicit `--threads` flag.

- [x] **`crossbeam-channel`** — bounded MPMC channels for pipeline stages.
  - Reference impl: — (Rust-native; analogous to Go channels)
  - Existing Rust: [`crossbeam-channel`](https://github.com/crossbeam-rs/crossbeam) `0.5.15`
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: Go channels; C++ `tbb::concurrent_bounded_queue`
  - Parallelism: this **is** the channel layer that lets stages run in parallel
  - SIMD: none
  - GPU-amenable: no
  - Upstream license: `MIT OR Apache-2.0`
  - Priority: `P0`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Quadrant ①. Adopt for reader → worker → writer pipelines (the standard bioinformatics streaming pattern). Avoid `std::sync::mpsc` (slower, SPSC only).

- [~] **`tokio` async I/O** — for cloud/network-bound workloads.
  - Reference impl: — (Rust-native; analogous to Go runtime)
  - Existing Rust: [`tokio`](https://github.com/tokio-rs/tokio) `1.52.3`
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: Boost.Asio (C++), Go runtime
  - Parallelism: M:N async runtime; multi-threaded executor
  - SIMD: none (runtime layer)
  - GPU-amenable: no
  - Upstream license: `MIT`
  - Priority: `P1`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Quadrant ①. noodles ships async variants of most readers behind a feature flag. Use sparingly — only when the workload is genuinely I/O-bound (htsget streaming, S3 / GCS BAM access). CPU-bound code does *not* benefit from async.

- [~] **`candle`** — pure-Rust deep-learning framework.
  - Reference impl: `Python/C++` · [pytorch/pytorch](https://github.com/pytorch/pytorch) · `BSD-3-Clause`
  - Existing Rust: [`candle-core`](https://crates.io/crates/candle-core) `0.10.2` (CPU, CUDA, Metal backends)
  - Existing Rust kind: `partial-port` (covers fewer ops than PyTorch but production-grade)
  - Existing non-C alternatives: PyTorch, JAX, MLX (Swift)
  - Parallelism: rayon for CPU tensor ops; CUDA / Metal for GPU
  - SIMD: explicit via the CPU backend's BLAS link (matmul); per-op kernels auto-vectorize
  - GPU-amenable: yes — primary purpose, with CUDA + Metal backends
  - Upstream license: `MIT OR Apache-2.0` (candle)
  - Priority: `P1`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Quadrant ①+②. CPU layer is pure-Rust + rayon; CUDA backend wraps `cudarc`. First-choice for DL inference (DeepVariant-class callers, AlphaFold). Smaller and faster-loading than PyTorch but covers fewer ops. Track `candle-transformers` for ESM/AlphaFold-relevant layers.

- [~] **`burn`** — alternative pure-Rust DL framework with backend swapping.
  - Reference impl: `Python/C++` · [pytorch/pytorch](https://github.com/pytorch/pytorch) · `BSD-3-Clause`
  - Existing Rust: [`burn`](https://github.com/tracel-ai/burn) `0.21.0` (`burn-wgpu`, `burn-candle`, `burn-ndarray`, `burn-tch` backends)
  - Existing Rust kind: `pure-port` (Rust orchestration; backend-dependent for kernels)
  - Existing non-C alternatives: PyTorch, TensorFlow
  - Parallelism: rayon CPU + wgpu / candle / tch backends
  - SIMD: explicit via backend (wgpu is GPU, ndarray is CPU SIMD via BLAS)
  - GPU-amenable: yes — cross-vendor via `burn-wgpu` (NVIDIA / AMD / Apple / WebGPU)
  - Upstream license: `MIT OR Apache-2.0`
  - Priority: `P1`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Quadrant ① (with wgpu / ndarray backends) or ② (with tch backend that wraps libtorch). More framework-y than candle (full training loop, autodiff, dashboard). The wgpu backend means cross-vendor GPU without a CUDA dependency — useful for portable inference binaries.

- [~] **`wgpu`** — portable GPU compute via WebGPU shading language.
  - Reference impl: — (browser standard; native via Vulkan/Metal/D3D12)
  - Existing Rust: [`wgpu`](https://github.com/gfx-rs/wgpu) `29.0.3`
  - Existing Rust kind: `pure-port` (wraps the OS graphics APIs, which is intrinsic to the layer)
  - Existing non-C alternatives: CUDA (NVIDIA only), Vulkan compute
  - Parallelism: GPU dispatch; CPU side is rayon-able
  - SIMD: GPU SIMT
  - GPU-amenable: yes — primary purpose
  - Upstream license: `MIT OR Apache-2.0`
  - Priority: `P2`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Quadrant ① at the orchestration layer (pure-Rust over Vulkan/Metal/D3D12 system APIs is by definition what a portable GPU library does). Low-level GPU access for custom kernels (Smith-Waterman, rolling hash) when DL frameworks are overkill. Adopt indirectly via `burn-wgpu`; direct use is a later-phase concern.

- [ ] **CUDA via `cust` / `cudarc`** — direct NVIDIA GPU compute.
  - Reference impl: `C++/CUDA` · [NVIDIA/cccl](https://github.com/NVIDIA/cccl) · `Apache-2.0`
  - Existing Rust: [`cust`](https://github.com/Rust-GPU/Rust-CUDA) `0.3.2`; [`cudarc`](https://github.com/coreylowman/cudarc) `0.19.6`
  - Existing Rust kind: `FFI-wrapper`
  - Existing non-C alternatives: CUDA C++, JAX
  - Parallelism: GPU SIMT
  - SIMD: GPU SIMT
  - GPU-amenable: yes — NVIDIA-only
  - Upstream license: `Apache-2.0` (cccl); `MIT OR Apache-2.0` for the Rust crates
  - Priority: `P2`
  - Layer: `adopt`
  - Consumes primitives: —
  - Notes: Quadrant ②. Only worth direct use when a tool ships a CUDA inner kernel (e.g. GPU BWA, GPU minimap2 forks). For DL inference, candle/burn already wrap `cudarc`. `cudarc` (coreylowman) is the lighter and more actively-maintained option; `cust` is broader but the Rust-GPU project moves slower.

- [ ] **SIMD portable abstractions** — `std::simd` and the `core_arch` intrinsics.
  - Reference impl: `C++` · intrinsics + `std::experimental::simd`
  - Existing Rust: `std::simd` (nightly, stabilising); [`wide`](https://github.com/Lokathor/wide) `1.4.0`; [`pulp`](https://github.com/sarah-quinones/pulp) `0.22.2` (runtime dispatch)
  - Existing Rust kind: `pure-port`
  - Existing non-C alternatives: ISPC; Highway (C++)
  - Parallelism: per-lane SIMD; combined with rayon for outer loop
  - SIMD: explicit — this **is** the SIMD layer
  - GPU-amenable: no — CPU SIMD only
  - Upstream license: `MIT OR Apache-2.0` (wide, pulp)
  - Priority: `P0`
  - Layer: `A` (foundation — every Layer A perf-critical crate consumes this)
  - Consumes primitives: —
  - Notes: Quadrant ①. Every hot kernel ships a scalar fallback and a SIMD path. `pulp` handles runtime feature detection nicely (`x86_64-v3` / `aarch64+neon` etc.). Avoid AVX-512-only code paths in shipped binaries (still poorly supported on consumer CPUs); detect at runtime via `pulp` instead.
