#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::missing_errors_doc,
    clippy::similar_names,
    clippy::many_single_char_names,
    clippy::float_cmp
)]

pub mod diversity;
pub mod hw;
pub mod ld;

pub use diversity::{tajimas_d, theta_pi, watterson_theta};
pub use hw::hwe_exact;
pub use ld::r_squared;

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum PopgenError {
    #[error("empty input")]
    Empty,
    #[error("sample size too small (n={n}, need ≥ {required})")]
    SampleTooSmall { n: usize, required: usize },
    #[error("invalid allele count: observed={observed}, expected ≤ {max}")]
    InvalidAlleleCount { observed: u64, max: u64 },
    #[error("no segregating sites — statistic undefined")]
    NoSegregating,
}

pub type Result<T> = std::result::Result<T, PopgenError>;
