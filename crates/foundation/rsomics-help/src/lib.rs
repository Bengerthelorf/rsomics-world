//! Family-wide `--help` renderer for `rsomics-*` CLIs.
//!
//! ## Usage
//!
//! Each binary declares a [`HelpSpec`] (typically as a `const`) and
//! calls [`render`] once from `main`. There is no per-binary rendering
//! code, no `println!` plumbing, no manual ANSI / layout. Rich / plain /
//! JSON modes are all dispatched from the same data.
//!
//! ```ignore
//! const HELP: rsomics_help::HelpSpec = rsomics_help::HelpSpec {
//!     name: env!("CARGO_PKG_NAME"),
//!     version: env!("CARGO_PKG_VERSION"),
//!     tagline: "what this tool does",
//!     origin: None,
//!     usage_lines: &["[OPTIONS] <INPUT>"],
//!     sections: &[rsomics_help::Section {
//!         title: "OPTIONS",
//!         flags: &[ /* FlagSpec ... */ ],
//!     }],
//!     examples: &[],
//!     json_result_schema_doc: None,
//! };
//!
//! fn main() {
//!     let args: Vec<String> = std::env::args().collect();
//!     if let Some(mode) = rsomics_help::intercept_help(&args) {
//!         rsomics_help::render(&HELP, mode);
//!         return;
//!     }
//!     // normal pipeline
//! }
//! ```
//!
//! ## Rendering contract
//!
//! - **Rich** (default on TTY): figlet banner + teal→violet gradient,
//!   bold colored section headers, green flag names, dim values /
//!   descriptions, two-line examples (`# description` then command).
//! - **Plain** (`--help --plain`, or auto on pipe / `NO_COLOR`): no
//!   banner, em-dash title line, ALL-CAPS section headers (no colon),
//!   no ANSI, one-line examples (`command   # description`).
//! - **JSON** (`--help --json`): structured envelope with
//!   `schema_version`, `tool`, `tool_version`, `tagline`, `origin`,
//!   `sections[].flags[]`, `examples`, `json_result_schema_doc`. The
//!   schema-version is owned by this crate, not the binary.

pub mod ansi;
pub mod banner;
pub mod modes;
pub mod render;
pub mod spec;

pub use ansi::no_color_env;
pub use banner::{Banner, Gradient};
pub use modes::{HelpMode, detect_mode, intercept_help};
pub use render::render;
pub use spec::{Example, FlagSpec, HELP_SCHEMA_VERSION, HelpSpec, Origin, Section};
