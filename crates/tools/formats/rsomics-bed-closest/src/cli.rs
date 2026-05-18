use clap::Parser;
use rsomics_bed_closest::closest;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-closest", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(short = 'a')]
    file_a: PathBuf,
    #[arg(short = 'b')]
    file_b: PathBuf,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        closest(&self.file_a, &self.file_b, &mut out)?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Find closest feature in B for each feature in A.",
    origin: Some(Origin {
        upstream: "bedtools closest",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["-a A.bed -b B.bed"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('a'),
                long: "a",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Input BED A.",
                why_default: None,
            },
            FlagSpec {
                short: Some('b'),
                long: "b",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Input BED B.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Closest gene to each peak",
        command: "rsomics-bed-closest -a peaks.bed -b genes.bed",
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
