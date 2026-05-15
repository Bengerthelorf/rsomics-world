use std::process;

use serde::Serialize;

use crate::error::Result;
use crate::exit::ExitCode;
use crate::flags::CommonFlags;
use crate::json::{ToolMeta, emit_error, emit_ok};
use crate::log::StderrLog;

pub fn run<T, F>(common: &CommonFlags, meta: ToolMeta, body: F) -> process::ExitCode
where
    F: FnOnce() -> Result<T>,
    T: Serialize,
{
    let log = StderrLog::from_flags(common);

    if let Err(e) = common.install_rayon_pool() {
        if common.json {
            emit_error(&meta, &e);
        }
        log.error(format_args!("{e}"));
        return ExitCode::from(&e).into();
    }

    match body() {
        Ok(result) => {
            if common.json {
                emit_ok(&meta, &result);
            }
            ExitCode::Ok.into()
        }
        Err(e) => {
            if common.json {
                emit_error(&meta, &e);
            }
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

    const META: ToolMeta = ToolMeta {
        name: "rsomics-runner-test",
        version: "0.0.0",
    };

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
        let code = run(&common, META, || Ok::<_, RsomicsError>(()));
        let expected: process::ExitCode = ExitCode::Ok.into();
        assert_eq!(format!("{code:?}"), format!("{expected:?}"));
    }

    #[test]
    fn invalid_input_body_maps_to_exit_one() {
        let common = cli();
        let code = run(&common, META, || -> Result<()> {
            Err(RsomicsError::InvalidInput("bad".into()))
        });
        let expected: process::ExitCode = ExitCode::InvalidInput.into();
        assert_eq!(format!("{code:?}"), format!("{expected:?}"));
    }

    #[test]
    fn rayon_pool_double_install_is_tolerated() {
        let common = cli();
        let code = run(&common, META, || Ok::<_, RsomicsError>(()));
        let expected: process::ExitCode = ExitCode::Ok.into();
        assert_eq!(format!("{code:?}"), format!("{expected:?}"));
    }
}
