use std::process;

use crate::error::Result;
use crate::exit::ExitCode;
use crate::flags::CommonFlags;
use crate::log::StderrLog;

/// Standard tool entrypoint. Every `rsomics-*` binary's `main` calls this:
///
/// ```ignore
/// fn main() -> std::process::ExitCode {
///     let args = Cli::parse();
///     rsomics_common::run(&args.common, || pipeline(args))
/// }
/// ```
///
/// Responsibilities, in order:
///
/// 1. Install the global rayon pool sized to `--threads`. A failure here
///    becomes a `ConfigError` and the process exits with [`ExitCode::ConfigError`]
///    before `body` runs.
/// 2. Run `body`. Whatever `Result` it returns is mapped to an [`ExitCode`].
/// 3. On error, print the error to stderr (unconditionally — errors are
///    not silenced by `--quiet`).
///
/// The `process::ExitCode` returned is what `main` should return directly.
pub fn run<F>(common: &CommonFlags, body: F) -> process::ExitCode
where
    F: FnOnce() -> Result<()>,
{
    let log = StderrLog::from_flags(common);

    if let Err(e) = common.install_rayon_pool() {
        log.error(format_args!("{e}"));
        return ExitCode::from(&e).into();
    }

    match body() {
        Ok(()) => ExitCode::Ok.into(),
        Err(e) => {
            log.error(format_args!("{e}"));
            ExitCode::from(&e).into()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::RsomicsError;
    use clap::Parser;

    #[derive(Parser)]
    struct Cli {
        #[command(flatten)]
        common: CommonFlags,
    }

    fn cli() -> CommonFlags {
        Cli::parse_from(["test", "--threads", "2"]).common
    }

    #[test]
    fn ok_body_exits_zero() {
        let common = cli();
        let code = run(&common, || Ok(()));
        // process::ExitCode doesn't expose its u8 directly; round-trip
        // through the public API to compare.
        let expected: process::ExitCode = ExitCode::Ok.into();
        assert_eq!(format!("{code:?}"), format!("{expected:?}"));
    }

    #[test]
    fn invalid_input_body_maps_to_exit_one() {
        let common = cli();
        let code = run(&common, || Err(RsomicsError::InvalidInput("bad".into())));
        let expected: process::ExitCode = ExitCode::InvalidInput.into();
        assert_eq!(format!("{code:?}"), format!("{expected:?}"));
    }

    #[test]
    fn rayon_pool_double_install_is_tolerated() {
        // `cli()` and `ok_body_exits_zero` will have already installed once
        // on a prior test; calling again must not surface a ConfigError.
        let common = cli();
        let code = run(&common, || Ok(()));
        let expected: process::ExitCode = ExitCode::Ok.into();
        assert_eq!(format!("{code:?}"), format!("{expected:?}"));
    }
}
