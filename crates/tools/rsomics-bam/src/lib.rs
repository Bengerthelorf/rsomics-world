//! Multi-subcommand `rsomics-*` template — each subcommand lives in its
//! own `src/cmd/<name>.rs` module and is dispatched from `main.rs`.
//! Appropriate when the tool has multiple distinct operating modes
//! (`view`, `sort`, `index`, …).
//!
//! For single-pipeline tools see [`rsomics-fastp`](https://crates.io/crates/rsomics-fastp)
//! — that crate's flat `Args` struct in `src/main.rs` is the right
//! template to clone when the tool does exactly one thing.

pub mod cmd;
pub mod htslib_bridge;
