use std::process;

use crate::error::RsomicsError;

/// Exit codes are part of the tool contract — callers branch on them.
/// New variants are additions, not renumberings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum ExitCode {
    Ok = 0,
    InvalidInput = 1,
    ConfigError = 2,
    UpstreamError = 3,
    IoError = 4,
}

impl From<&RsomicsError> for ExitCode {
    fn from(e: &RsomicsError) -> Self {
        match e {
            RsomicsError::Io(_) => Self::IoError,
            RsomicsError::InvalidInput(_) => Self::InvalidInput,
            RsomicsError::ConfigError(_) => Self::ConfigError,
            RsomicsError::UpstreamError(_) => Self::UpstreamError,
        }
    }
}

impl From<ExitCode> for process::ExitCode {
    fn from(code: ExitCode) -> Self {
        Self::from(code as u8)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ok_is_zero() {
        assert_eq!(ExitCode::Ok as u8, 0);
    }

    #[test]
    fn each_error_variant_maps_to_distinct_code() {
        let io = RsomicsError::Io(std::io::Error::other("x"));
        let inv = RsomicsError::InvalidInput("x".into());
        let cfg = RsomicsError::ConfigError("x".into());
        let up = RsomicsError::UpstreamError("x".into());
        assert_eq!(ExitCode::from(&io), ExitCode::IoError);
        assert_eq!(ExitCode::from(&inv), ExitCode::InvalidInput);
        assert_eq!(ExitCode::from(&cfg), ExitCode::ConfigError);
        assert_eq!(ExitCode::from(&up), ExitCode::UpstreamError);
    }

    #[test]
    fn converts_into_process_exit_code() {
        let _: process::ExitCode = ExitCode::InvalidInput.into();
    }
}
