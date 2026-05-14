//! Section / row primitives shared across binaries.
//!
//! These return owned `String`s so the binary's `--help` body composes
//! naturally with `println!`. None of them allocate global state; each
//! call is independent and `color` controls whether ANSI escapes are
//! emitted.

use crate::ansi::{Palette, bold, dim, rgb};

/// Bold section header like `USAGE:` / `OPTIONS:` / `EXAMPLES:`.
#[must_use]
pub fn section_header(title: &str, color: bool) -> String {
    bold(color, &format!("{title}:"))
}

/// Tagline shown directly under the banner:
/// `name v0.X  ─  description`.
#[must_use]
pub fn tagline(name: &str, version: &str, description: &str, color: bool) -> String {
    let name_part = bold(color, &rgb(color, Palette::TEAL, name));
    let version_part = dim(color, &format!("v{version}"));
    let sep = dim(color, "─");
    let desc_part = dim(color, description);
    format!("{name_part} {version_part}  {sep}  {desc_part}")
}

/// One flag-table row. `value` is the metavariable (e.g. `<path>`) and
/// `desc` is the help text. Set `short` to `None` when the flag has
/// no short alias.
///
/// `desc_col` is the absolute column at which descriptions start. The
/// row pads (or, if the flag+value is already too wide, just inserts
/// two spaces) so descriptions line up in a column regardless of how
/// long `--long_name` or `<value>` happens to be.
#[must_use]
pub fn flag_row(
    short: Option<char>,
    long: &str,
    value: Option<&str>,
    desc: &str,
    color: bool,
    desc_col: usize,
) -> String {
    // Visible-width budget (ignoring ANSI escapes since they don't
    // contribute to terminal columns):
    //   2 indent + short_part (4 if short, 4 spaces otherwise)
    //   + "--" + long
    //   + (" " + value) if value
    let short_visible = match short {
        Some(_) => "-x, ".len(),
        None => "    ".len(),
    };
    let value_visible = value.map_or(0, |v| 1 + v.len());
    let used = 2 + short_visible + 2 + long.len() + value_visible;
    let pad = desc_col.saturating_sub(used).max(2);

    let short_part = match short {
        Some(c) => format!("{}, ", rgb(color, Palette::GREEN, &format!("-{c}"))),
        None => "    ".to_string(),
    };
    let long_painted = rgb(color, Palette::GREEN, &format!("--{long}"));
    let value_part = value.map_or_else(String::new, |v| format!(" {}", dim(color, v)));
    let desc_part = dim(color, desc);
    format!(
        "  {short_part}{long_painted}{value_part}{}{desc_part}",
        " ".repeat(pad)
    )
}

/// Per-row spec for [`flag_table`].
pub struct FlagRowSpec<'a> {
    pub short: Option<char>,
    pub long: &'a str,
    pub value: Option<&'a str>,
    pub desc: &'a str,
}

/// Render a group of flags with descriptions aligned at a common column
/// computed from the widest flag in the group.
#[must_use]
pub fn flag_table(rows: &[FlagRowSpec<'_>], color: bool) -> String {
    let widest = rows
        .iter()
        .map(|r| {
            // 2 indent + "-x, " (4) or "    " (4) + "--" (2) + long + value
            2 + 4 + 2 + r.long.len() + r.value.map_or(0, |v| 1 + v.len())
        })
        .max()
        .unwrap_or(0);
    let desc_col = widest + 2;
    rows.iter()
        .map(|r| flag_row(r.short, r.long, r.value, r.desc, color, desc_col))
        .collect::<Vec<_>>()
        .join("\n")
}

/// Two-line example: dim description followed by indented command.
#[must_use]
pub fn example_line(description: &str, command: &str, color: bool) -> String {
    let desc = dim(color, &format!("  # {description}"));
    let cmd = rgb(color, Palette::GREEN, &format!("  {command}"));
    format!("{desc}\n{cmd}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn section_header_bolds_in_color_mode() {
        let s = section_header("USAGE", true);
        assert!(s.contains("USAGE:"));
        assert!(s.contains("\x1b["));
    }

    #[test]
    fn section_header_plain_has_no_escape() {
        let s = section_header("USAGE", false);
        assert_eq!(s, "USAGE:");
    }

    #[test]
    fn flag_row_with_short_includes_short() {
        let s = flag_row(Some('i'), "in1", Some("<path>"), "R1 input", false, 16);
        assert!(s.contains("-i"));
        assert!(s.contains("--in1"));
        assert!(s.contains("<path>"));
        assert!(s.contains("R1 input"));
    }

    #[test]
    fn flag_row_without_short_pads_correctly() {
        let s = flag_row(
            None,
            "trim_tail1",
            None,
            "Bases trimmed from R1 3'",
            false,
            16,
        );
        assert!(s.contains("--trim_tail1"));
        assert!(s.contains("Bases trimmed"));
    }
}
