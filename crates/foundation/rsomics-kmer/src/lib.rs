#![allow(
    clippy::cast_possible_truncation,
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::must_use_candidate
)]

pub mod count;
pub mod encode;
pub mod hash;
pub mod iter;

pub use count::KmerCounts;
pub use encode::{Kmer, canonical, decode, encode, reverse_complement};
pub use hash::{murmur3_x64_128, nthash_iter, nthash_one};
pub use iter::KmerIter;

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum KmerError {
    #[error("k must be in 1..=32 (got {0})")]
    KOutOfRange(usize),
    #[error("sequence shorter than k: len={len}, k={k}")]
    SeqTooShort { len: usize, k: usize },
    #[error("non-ACGT base at position {pos}: {byte:?}")]
    NonAcgt { pos: usize, byte: u8 },
}

pub type Result<T> = std::result::Result<T, KmerError>;
