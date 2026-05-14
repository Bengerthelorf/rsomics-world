use std::io;

use thiserror::Error;

/// Canonical error type for every `rsomics-*` tool. Variants carve out the
/// distinctions that callers actually act on: input validity, configuration
/// validity, upstream-tool / FFI failures, and I/O. Anything else collapses
/// into one of these; we deliberately keep the variant count small so the
/// [`crate::ExitCode`] mapping stays exhaustive without ceremony.
#[derive(Debug, Error)]
pub enum RsomicsError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),

    #[error("invalid input: {0}")]
    InvalidInput(String),

    #[error("configuration error: {0}")]
    ConfigError(String),

    #[error("upstream tool error: {0}")]
    UpstreamError(String),
}

/// Shorthand `Result` parameterised over [`RsomicsError`]. Tool code should
/// return this from any fallible boundary; the runner ([`crate::run`]) maps
/// it to a process exit code.
pub type Result<T> = std::result::Result<T, RsomicsError>;

impl From<std::num::ParseIntError> for RsomicsError {
    fn from(e: std::num::ParseIntError) -> Self {
        Self::InvalidInput(e.to_string())
    }
}

impl From<std::num::ParseFloatError> for RsomicsError {
    fn from(e: std::num::ParseFloatError) -> Self {
        Self::InvalidInput(e.to_string())
    }
}

impl From<std::str::Utf8Error> for RsomicsError {
    fn from(e: std::str::Utf8Error) -> Self {
        Self::InvalidInput(e.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_includes_message() {
        let e = RsomicsError::InvalidInput("bad header".into());
        assert_eq!(e.to_string(), "invalid input: bad header");
    }

    #[test]
    fn io_from_conversion_round_trips() {
        let underlying = io::Error::new(io::ErrorKind::NotFound, "missing.fastq");
        let e: RsomicsError = underlying.into();
        match e {
            RsomicsError::Io(_) => {}
            other => panic!("expected Io variant, got {other:?}"),
        }
    }

    #[test]
    fn parse_int_routes_to_invalid_input() {
        let e: RsomicsError = "abc".parse::<u32>().unwrap_err().into();
        assert!(matches!(e, RsomicsError::InvalidInput(_)));
    }
}
