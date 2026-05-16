#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::many_single_char_names,
    clippy::float_cmp,
    clippy::similar_names
)]

pub mod combine;
pub mod fdr;
pub mod hypothesis;

pub use combine::{fisher_combine, stouffer_combine};
pub use fdr::{
    bh_adjust, bonferroni_adjust, by_adjust, hochberg_adjust, holm_adjust, hommel_adjust,
    none_adjust,
};
pub use hypothesis::{Alternative, fisher_exact_2x2, mann_whitney_u, welch_t};

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum StatsError {
    #[error("empty input")]
    Empty,
    #[error("sample too small (n={n}, need ≥ {required})")]
    SampleTooSmall { n: usize, required: usize },
    #[error("zero variance in both samples — t statistic is undefined")]
    ZeroVariance,
    #[error("p-value out of range: {0}")]
    InvalidPValue(f64),
    #[error("statrs error: {0}")]
    Statrs(String),
}

pub type Result<T> = std::result::Result<T, StatsError>;
