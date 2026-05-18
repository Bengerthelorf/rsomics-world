use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_gff_filter::filter_gff;
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-gff-filter", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(value_name = "INPUT")]
    input: PathBuf,
    /// Keep only features of this type (column 3: gene, exon, CDS, etc.).
    #[arg(long = "type")]
    feature_type: Option<String>,
    /// Keep lines matching this regex (searched across the full line).
    #[arg(short = 'p', long = "pattern")]
    pattern: Option<String>,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        filter_gff(
            &self.input,
            self.feature_type.as_deref(),
            self.pattern.as_deref(),
            &mut out,
        )?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Filter GFF/GTF annotations by type or pattern.",
    origin: Some(Origin {
        upstream: "awk on GFF / AGAT",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["[--type TYPE] [-p PATTERN] <INPUT.gff>"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "INPUT",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Input GFF/GTF.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "type",
                aliases: &[],
                value: Some("<str>"),
                type_hint: Some("String"),
                required: false,
                default: None,
                description: "Feature type (gene, exon, CDS...).",
                why_default: None,
            },
            FlagSpec {
                short: Some('p'),
                long: "pattern",
                aliases: &[],
                value: Some("<regex>"),
                type_hint: Some("String"),
                required: false,
                default: None,
                description: "Regex pattern.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Extract only exons",
        command: "rsomics-gff-filter --type exon genes.gff > exons.gff",
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
