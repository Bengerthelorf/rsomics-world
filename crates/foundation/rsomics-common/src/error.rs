use std::io;

use thiserror::Error;

#[derive(Debug, Error)]
#[non_exhaustive]
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

/// Attach a contextual prefix to a `Result`.
pub trait Context<T> {
    #[allow(clippy::missing_errors_doc)]
    fn rs_context(self, msg: impl Into<String>) -> Result<T>;

    #[allow(clippy::missing_errors_doc)]
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

impl<T> Context<T> for Result<T> {
    fn rs_context(self, msg: impl Into<String>) -> Result<T> {
        self.map_err(|e| prepend(&msg.into(), e))
    }

    fn rs_with_context<F>(self, f: F) -> Result<T>
    where
        F: FnOnce() -> String,
    {
        self.map_err(|e| prepend(&f(), e))
    }
}

fn prepend(prefix: &str, e: RsomicsError) -> RsomicsError {
    match e {
        RsomicsError::Io(inner) => {
            let kind = inner.kind();
            RsomicsError::Io(io::Error::new(kind, format!("{prefix}: {inner}")))
        }
        RsomicsError::InvalidInput(s) => RsomicsError::InvalidInput(format!("{prefix}: {s}")),
        RsomicsError::ConfigError(s) => RsomicsError::ConfigError(format!("{prefix}: {s}")),
        RsomicsError::UpstreamError(s) => RsomicsError::UpstreamError(format!("{prefix}: {s}")),
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

    #[test]
    fn rs_context_chains_through_rsomics_error_and_preserves_variant() {
        let res: Result<()> = Err(RsomicsError::InvalidInput("bad header".into()));
        let err = res.rs_context("parsing record 17").unwrap_err();
        assert_eq!(
            err.to_string(),
            "invalid input: parsing record 17: bad header"
        );
        assert!(matches!(err, RsomicsError::InvalidInput(_)));
    }

    #[test]
    fn rs_context_chains_through_io_variant_preserves_kind() {
        let inner = io::Error::new(io::ErrorKind::PermissionDenied, "denied");
        let res: Result<()> = Err(RsomicsError::Io(inner));
        let err = res.rs_context("opening output").unwrap_err();
        let RsomicsError::Io(io_err) = err else {
            panic!("expected Io variant");
        };
        assert_eq!(io_err.kind(), io::ErrorKind::PermissionDenied);
    }
}
