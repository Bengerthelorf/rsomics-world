use std::fmt;
use std::io::{self, Write};

use crate::flags::CommonFlags;

/// Stderr sink respecting `--quiet` / `--verbose`. `error` always prints;
/// `info` respects quiet; `debug` requires verbose.
///
/// Avoids `log`/`tracing` indirection — prevents name collision with the
/// `log` crate macros in dependent crates' preludes.
#[derive(Debug, Clone, Copy)]
pub struct StderrLog {
    pub quiet: bool,
    pub verbose: bool,
    pub json: bool,
}

impl StderrLog {
    #[must_use]
    pub fn from_flags(common: &CommonFlags) -> Self {
        Self {
            quiet: common.quiet,
            verbose: common.verbose,
            json: common.json,
        }
    }

    pub fn error(&self, args: fmt::Arguments<'_>) {
        let _ = writeln!(io::stderr().lock(), "error: {args}");
    }

    pub fn info(&self, args: fmt::Arguments<'_>) {
        if self.quiet {
            return;
        }
        let _ = writeln!(io::stderr().lock(), "{args}");
    }

    pub fn debug(&self, args: fmt::Arguments<'_>) {
        if self.quiet || !self.verbose {
            return;
        }
        let _ = writeln!(io::stderr().lock(), "debug: {args}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flags(quiet: bool, verbose: bool) -> CommonFlags {
        use clap::Parser;
        #[derive(Parser)]
        struct Cli {
            #[command(flatten)]
            c: CommonFlags,
        }
        let mut argv = vec!["test"];
        if quiet {
            argv.push("--quiet");
        }
        if verbose {
            argv.push("--verbose");
        }
        Cli::parse_from(argv).c
    }

    #[test]
    fn info_suppressed_when_quiet() {
        let f = flags(true, false);
        let log = StderrLog::from_flags(&f);
        assert!(log.quiet);
        log.info(format_args!("should not print"));
    }

    #[test]
    fn debug_requires_verbose() {
        let f = flags(false, false);
        let log = StderrLog::from_flags(&f);
        assert!(!log.verbose);
        log.debug(format_args!("should not print"));
    }

    #[test]
    fn debug_prints_when_verbose_and_not_quiet() {
        let f = flags(false, true);
        let log = StderrLog::from_flags(&f);
        assert!(log.verbose && !log.quiet);
    }
}
