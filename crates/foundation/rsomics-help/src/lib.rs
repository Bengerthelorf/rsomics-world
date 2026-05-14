//! Family-wide `--help` renderer for `rsomics-*` CLIs.
//!
//! - Banner: figlet-rs auto-generates the ASCII art from the binary
//!   name; this crate adds a per-line gradient.
//! - Mode detection: `--plain` / `--json` / `NO_COLOR` env / stdout-isn't-a-tty
//!   each fold into a [`HelpMode`].
//! - Section / flag helpers keep the rendering across ~150 future
//!   binaries visually uniform without forcing a rigid schema.
//!
//! Layout philosophy: this crate provides primitives, not a framework.
//! Each binary still composes its own `--help` body. The crate ensures
//! the *look* is consistent; the *content* stays in the binary.

pub mod ansi;
pub mod banner;
pub mod json;
pub mod modes;
pub mod section;

pub use ansi::{Palette, no_color_env};
pub use banner::{Banner, Gradient};
pub use json::{Example, FlagGroup, FlagSpec, HelpJson, Origin};
pub use modes::{HelpMode, detect_mode, intercept_help};
pub use section::{FlagRowSpec, example_line, flag_row, flag_table, section_header, tagline};
