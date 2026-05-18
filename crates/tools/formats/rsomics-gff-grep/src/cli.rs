use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_gff_grep::grep;
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-gff-grep", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(short = 'p', long = "pattern")]
    pattern: String,
    #[arg(value_name = "INPUT")]
    input: PathBuf,
    /// Search only the attributes column (col 9).
    #[arg(long = "attr-only")]
    attr_only: bool,
    #[arg(long = "invert-match")]
    invert: bool,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        grep(
            &self.input,
            &self.pattern,
            self.attr_only,
            self.invert,
            &mut out,
        )?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Filter GFF/GTF by attribute regex.",
    origin: Some(Origin {
        upstream: "grep / AGAT",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["-p PATTERN [--attr-only] <INPUT.gff>"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('p'),
                long: "pattern",
                aliases: &[],
                value: Some("<regex>"),
                type_hint: Some("String"),
                required: true,
                default: None,
                description: "Regex pattern.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "attr-only",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Search attributes column only.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "invert-match",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Output non-matching.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Extract TP53 gene features",
        command: "rsomics-gff-grep -p 'gene_name=TP53' --attr-only genes.gff",
    }],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;
    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
