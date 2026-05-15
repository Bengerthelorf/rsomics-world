use std::io::IsTerminal;

use crate::ansi::no_color_env;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum HelpMode {
    Rich,
    Plain,
    Json,
}

#[must_use]
pub fn wants_help(args: &[String]) -> bool {
    args.iter()
        .skip(1)
        .any(|a| a == "--help" || a == "-h" || a == "help")
}

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
