#![allow(clippy::missing_errors_doc)]

pub mod dedup;
pub mod report;

pub use dedup::{DedupConfig, DedupMode, run_se};
pub use report::DedupReport;
