use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_fasta_grep::grep;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fasta-grep", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(short = 'p', long = "pattern")]
    pattern: String,

    #[arg(value_name = "INPUT")]
    input: PathBuf,

    #[arg(long = "invert-match")]
    invert: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        grep(&self.input, &self.pattern, self.invert, &mut out)?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Filter FASTA records by name regex.",
    origin: Some(Origin {
        upstream: "seqkit grep (FASTA mode)",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["-p PATTERN <INPUT.fa>"],
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
                description: "Regex to match sequence names.",
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
                description: "Output non-matching sequences.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Extract chromosome 1",
        command: "rsomics-fasta-grep -p '^chr1$' genome.fa > chr1.fa",
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
