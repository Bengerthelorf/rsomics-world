use std::io;

use thiserror::Error;

/// Canonical error type for every `rsomics-*` tool. Variants carve out the
/// distinctions that callers actually act on: input validity, configuration
/// validity, upstream-tool / FFI failures, and I/O. Anything else collapses
/// into one of these; we deliberately keep the variant count small so the
/// [`crate::ExitCode`] mapping stays exhaustive without ceremony.
#[derive(Debug, Error)]
pub enum RsomicsError {
    /// Any failure surfaced through `std::io::Error`. Note this also catches
    /// "user gave a non-existent path" cases — at this layer we can't tell
    /// a hardware fault from a typo, so both route here. Integrators that
    /// need finer granularity should bracket their boundary code with an
    /// explicit `Path::exists()` check and map to `InvalidInput`.
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

/// Attach a description to a `Result` that turns into the prefix of the
/// resulting `RsomicsError` message. Lets tools write
/// `File::create(p).rs_context(|| format!("opening {}", p.display()))?` and
/// get a contextual error without an `anyhow` dependency.
pub trait Context<T> {
    /// Eager-evaluated context. Use when the message is a literal or
    /// already-computed string.
    ///
    /// # Errors
    ///
    /// Returns the original `Err` wrapped as [`RsomicsError`] with the
    /// supplied prefix prepended.
    fn rs_context(self, msg: impl Into<String>) -> Result<T>;

    /// Lazy-evaluated context. The closure runs only when the result is
    /// `Err`, so `format!` work is avoided on the success path.
    ///
    /// # Errors
    ///
    /// Same as [`Self::rs_context`].
    fn rs_with_context<F>(self, f: F) -> Result<T>
    where
        F: FnOnce() -> String;
}

impl<T> Context<T> for std::result::Result<T, io::Error> {
    fn rs_context(self, msg: impl Into<String>) -> Result<T> {
        self.map_err(|e| {
            let kind = e.kind();
            RsomicsError::Io(io::Error::new(kind, format!("{}: {e}", msg.into())))
        })
    }

    fn rs_with_context<F>(self, f: F) -> Result<T>
    where
        F: FnOnce() -> String,
    {
        self.map_err(|e| {
            let kind = e.kind();
            RsomicsError::Io(io::Error::new(kind, format!("{}: {e}", f())))
        })
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

    #[test]
    fn rs_context_prefixes_error_message() {
        let res: std::result::Result<(), io::Error> =
            Err(io::Error::new(io::ErrorKind::NotFound, "missing.fastq"));
        let err = res.rs_context("opening input").unwrap_err();
        assert_eq!(err.to_string(), "I/O error: opening input: missing.fastq");
    }

    #[test]
    fn rs_with_context_is_lazy_on_ok() {
        let mut called = false;
        let res: std::result::Result<u32, io::Error> = Ok(7);
        let _ = res.rs_with_context(|| {
            called = true;
            "should not be evaluated".into()
        });
        assert!(!called, "closure must not run on Ok");
    }
}
