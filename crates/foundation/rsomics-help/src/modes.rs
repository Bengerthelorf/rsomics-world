//! `--help` mode detection.
//!
//! Priority order:
//! 1. `--json` → [`HelpMode::Json`] (machine-readable, AI-consumer path)
//! 2. `--plain` → [`HelpMode::Plain`] (clap-style, no banner / no colors)
//! 3. `NO_COLOR` env or stdout-isn't-a-tty → [`HelpMode::Plain`]
//! 4. Default → [`HelpMode::Rich`]
//!
//! [`intercept_help`] folds the "is the user asking for help at all?"
//! check on top so a binary's `main` can do:
//!
//! ```ignore
//! if let Some(mode) = rsomics_help::intercept_help(&std::env::args().collect::<Vec<_>>()) {
//!     match mode {
//!         HelpMode::Rich => print_rich(),
//!         HelpMode::Plain => Cli::command().print_help().unwrap(),
//!         HelpMode::Json => print_json(),
//!     }
//!     return ExitCode::SUCCESS;
//! }
//! let args = Cli::parse();
//! ```

use std::io::IsTerminal;

use crate::ansi::no_color_env;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HelpMode {
    Rich,
    Plain,
    Json,
}

/// True if any of `--help`, `-h`, or bare `help` appears in `args`.
/// Skips `args[0]` (the binary name).
#[must_use]
pub fn wants_help(args: &[String]) -> bool {
    args.iter()
        .skip(1)
        .any(|a| a == "--help" || a == "-h" || a == "help")
}

/// Pick the help mode given the full argv and the runtime environment.
#[must_use]
pub fn detect_mode(args: &[String]) -> HelpMode {
    let json = args.iter().any(|a| a == "--json");
    let plain = args.iter().any(|a| a == "--plain");
    if json {
        return HelpMode::Json;
    }
    if plain || no_color_env() || !std::io::stdout().is_terminal() {
        return HelpMode::Plain;
    }
    HelpMode::Rich
}

/// Convenience wrapper: returns `Some(mode)` iff the user is asking
/// for help. Always pair `wants_help` and `detect_mode`.
#[must_use]
pub fn intercept_help(args: &[String]) -> Option<HelpMode> {
    if wants_help(args) {
        Some(detect_mode(args))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(a: &[&str]) -> Vec<String> {
        a.iter().map(|s| (*s).to_string()).collect()
    }

    #[test]
    fn wants_help_picks_up_help_long() {
        assert!(wants_help(&argv(&["bin", "--help"])));
    }

    #[test]
    fn wants_help_picks_up_help_short() {
        assert!(wants_help(&argv(&["bin", "-h"])));
    }

    #[test]
    fn wants_help_picks_up_help_subcommand() {
        assert!(wants_help(&argv(&["bin", "help"])));
    }

    #[test]
    fn no_help_when_no_flag() {
        assert!(!wants_help(&argv(&["bin", "-i", "x.fq"])));
    }

    #[test]
    fn json_beats_plain() {
        let mode = detect_mode(&argv(&["bin", "--help", "--json", "--plain"]));
        assert_eq!(mode, HelpMode::Json);
    }
}
