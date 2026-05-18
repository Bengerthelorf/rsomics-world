use clap::Parser;
use rsomics_bed_to_gff::bed_to_gff;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-to-gff", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(value_name = "INPUT")]
    input: PathBuf,
    #[arg(long = "source", default_value = "rsomics")]
    source: String,
    #[arg(long = "type", default_value = "region")]
    feature_type: String,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        bed_to_gff(&self.input, &self.source, &self.feature_type, &mut out)?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Convert BED intervals to GFF3.",
    origin: Some(Origin {
        upstream: "bedtools / bed2gff",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["<INPUT.bed>"],
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
                description: "Input BED.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "source",
                aliases: &[],
                value: Some("<str>"),
                type_hint: Some("String"),
                required: false,
                default: Some("rsomics"),
                description: "GFF source field.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "type",
                aliases: &[],
                value: Some("<str>"),
                type_hint: Some("String"),
                required: false,
                default: Some("region"),
                description: "GFF feature type.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Convert peaks to GFF",
        command: "rsomics-bed-to-gff peaks.bed > peaks.gff",
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
