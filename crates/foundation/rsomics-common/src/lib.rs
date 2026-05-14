//! Shared primitives for every `rsomics-*` binary.
//!
//! Tools in this workspace share a small set of concerns — error type,
//! exit codes, CLI flags, stderr noise control, and the main-function
//! plumbing that ties them together. Each gets one module here, re-exported
//! from the crate root so the import line in any `main.rs` stays short.
//!
//! ## What's here
//!
//! - [`RsomicsError`] / [`Result`] — canonical error type, four variants
//!   covering input validity, configuration, upstream/FFI failures, and I/O.
//! - [`ExitCode`] — numeric process exit codes; the `From<&RsomicsError>`
//!   impl is the single source of truth for "which error becomes which
//!   exit code."
//! - [`CommonFlags`] — `--threads / --json / --quiet / --verbose / --seed`,
//!   flattened into each tool's `clap` struct.
//! - [`StderrLog`] — minimal quiet/verbose-aware stderr sink. `error!`
//!   always prints; `info!` honours `--quiet`; `debug!` requires `--verbose`.
//! - [`run`] — the entrypoint helper. Every tool's `main` calls this.
//!
//! ## Usage shape
//!
//! ```ignore
//! use clap::Parser;
//! use rsomics_common::{CommonFlags, Result, run};
//!
//! #[derive(Parser)]
//! struct Cli {
//!     /// tool-specific flags here
//!     #[arg(short = 'i', long)]
//!     input: std::path::PathBuf,
//!
//!     #[command(flatten)]
//!     common: CommonFlags,
//! }
//!
//! fn pipeline(args: Cli) -> Result<()> {
//!     // ...real work...
//!     Ok(())
//! }
//!
//! fn main() -> std::process::ExitCode {
//!     let args = Cli::parse();
//!     let common = args.common.clone();
//!     run(&common, || pipeline(args))
//! }
//! ```
//!
//! That layout gives every tool: thread-count negotiation, deterministic
//! seeding, a JSON-output flag wired up, a quiet/verbose contract for
//! stderr, and an exit-code mapping consistent with every other tool in
//! the family.

pub mod error;
pub mod exit;
pub mod flags;
pub mod log;
pub mod runner;

pub use error::{Context, Result, RsomicsError};
pub use exit::ExitCode;
pub use flags::CommonFlags;
pub use log::StderrLog;
pub use runner::run;
