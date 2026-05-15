use serde::Serialize;

use crate::ansi::{Palette, bold, dim, no_color_env, rgb};
use crate::banner::Banner;
use crate::modes::HelpMode;
use crate::spec::{FlagSpec, HELP_SCHEMA_VERSION, HelpSpec};

pub fn render(spec: &HelpSpec, mode: HelpMode) {
    match mode {
        HelpMode::Rich => render_rich(spec),
        HelpMode::Plain => render_plain(spec),
        HelpMode::Json => render_json(spec),
    }
}

fn render_rich(spec: &HelpSpec) {
    let color = !no_color_env();
    let banner = Banner::family(spec.name).render(color);
    if !banner.is_empty() {
        println!();
        println!("{banner}");
    }
    println!();
    println!(
        "  {}",
        rich_tagline(spec.name, spec.version, spec.tagline, color)
    );
    println!();
    println!("{}", bold(color, "USAGE:"));
    for line in spec.usage_lines {
        println!("  {} {line}", rgb(color, Palette::GREEN, spec.name));
    }
    for section in spec.sections {
        println!();
        println!("{}", bold(color, &format!("{}:", section.title)));
        println!("{}", render_flag_table(section.flags, color));
    }
    if !spec.examples.is_empty() {
        println!();
        println!("{}", bold(color, "EXAMPLES:"));
        for ex in spec.examples {
            println!("  {}", dim(color, &format!("# {}", ex.description)));
            println!("  {}", rgb(color, Palette::GREEN, ex.command));
        }
    }
    println!();
}

fn render_plain(spec: &HelpSpec) {
    println!("{} {} — {}", spec.name, spec.version, spec.tagline);
    println!();
    println!("USAGE");
    for line in spec.usage_lines {
        println!("  {} {line}", spec.name);
    }
    for section in spec.sections {
        println!();
        println!("{}", section.title);
        println!("{}", render_flag_table(section.flags, false));
    }
    if !spec.examples.is_empty() {
        println!();
        println!("EXAMPLES");
        let widest_cmd = spec
            .examples
            .iter()
            .map(|e| e.command.len())
            .max()
            .unwrap_or(0);
        for ex in spec.examples {
            let pad = widest_cmd.saturating_sub(ex.command.len()) + 4;
            println!("  {}{}# {}", ex.command, " ".repeat(pad), ex.description);
        }
    }
    println!();
}

#[derive(Serialize)]
struct JsonEnvelope<'a> {
    schema_version: &'static str,
    #[serde(flatten)]
    spec: &'a HelpSpec<'a>,
}

fn render_json(spec: &HelpSpec) {
    let env = JsonEnvelope {
        schema_version: HELP_SCHEMA_VERSION,
        spec,
    };
    let _ = serde_json::to_writer_pretty(std::io::stdout().lock(), &env);
    println!();
}

fn rich_tagline(name: &str, version: &str, description: &str, color: bool) -> String {
    let name_part = bold(color, &rgb(color, Palette::TEAL, name));
    let version_part = dim(color, &format!("v{version}"));
    let sep = dim(color, "─");
    let desc_part = dim(color, description);
    format!("{name_part} {version_part}  {sep}  {desc_part}")
}

fn render_flag_table(flags: &[FlagSpec<'_>], color: bool) -> String {
    let widest_visible = flags.iter().map(visible_flag_width).max().unwrap_or(0);
    let desc_col = widest_visible + 2;
    flags
        .iter()
        .map(|f| render_flag_row(f, color, desc_col))
        .collect::<Vec<_>>()
        .join("\n")
}

fn visible_flag_width(f: &FlagSpec<'_>) -> usize {
    // 2 indent + "-x, " (4) or "    " (4) + "--" (2) + long + value
    2 + 4 + 2 + f.long.len() + f.value.map_or(0, |v| 1 + v.len())
}

fn render_flag_row(f: &FlagSpec<'_>, color: bool, desc_col: usize) -> String {
    let used = visible_flag_width(f);
    let pad = desc_col.saturating_sub(used).max(2);
    let short_part = match f.short {
        Some(c) => format!("{}, ", rgb(color, Palette::GREEN, &format!("-{c}"))),
        None => "    ".to_string(),
    };
    let long_painted = rgb(color, Palette::GREEN, &format!("--{}", f.long));
    let value_part = f
        .value
        .map_or_else(String::new, |v| format!(" {}", dim(color, v)));
    let desc_part = dim(color, f.description);
    format!(
        "  {short_part}{long_painted}{value_part}{}{desc_part}",
        " ".repeat(pad)
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::{Example, Section};

    const SAMPLE: HelpSpec = HelpSpec {
        name: "rsomics-test",
        version: "0.0.0",
        tagline: "test tool",
        origin: None,
        usage_lines: &["[OPTIONS] <FILE>"],
        sections: &[Section {
            title: "OPTIONS",
            flags: &[
                FlagSpec {
                    short: Some('a'),
                    long: "all",
                    aliases: &[],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "do everything",
                    why_default: None,
                },
                FlagSpec {
                    short: None,
                    long: "verbose-output",
                    aliases: &[],
                    value: Some("<n>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: Some("0"),
                    description: "verbosity",
                    why_default: None,
                },
            ],
        }],
        examples: &[Example {
            description: "basic",
            command: "rsomics-test input.fa",
        }],
        json_result_schema_doc: None,
    };

    #[test]
    fn json_envelope_has_schema_version() {
        let env = JsonEnvelope {
            schema_version: HELP_SCHEMA_VERSION,
            spec: &SAMPLE,
        };
        let s = serde_json::to_string(&env).unwrap();
        assert!(s.contains("\"schema_version\":\"1.0\""));
        assert!(s.contains("\"tool\""), "{s}");
    }

    #[test]
    fn plain_table_aligns_to_widest() {
        let s = render_flag_table(SAMPLE.sections[0].flags, false);
        // Two rows; both descriptions should start at the same column.
        let lines: Vec<&str> = s.lines().collect();
        let col = |line: &str, needle: &str| line.find(needle).unwrap();
        assert_eq!(col(lines[0], "do everything"), col(lines[1], "verbosity"));
    }
}
