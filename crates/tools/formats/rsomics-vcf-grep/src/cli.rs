use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use rsomics_vcf_grep::grep;
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-vcf-grep", version, about, long_about = None, disable_help_flag = true)]
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
    tagline: "Filter VCF variants by regex pattern.",
    origin: Some(Origin {
        upstream: "bcftools view -i (regex subset)",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["-p PATTERN <INPUT.vcf>"],
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
                long: "invert-match",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Output non-matching variants.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Filter by rsID",
        command: "rsomics-vcf-grep -p 'rs[0-9]+' variants.vcf > with_rsid.vcf",
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
