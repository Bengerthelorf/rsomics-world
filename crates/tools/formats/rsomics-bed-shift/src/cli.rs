use clap::Parser;
use rsomics_bed_shift::shift;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-shift", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(short = 'i', long = "input")]
    input: PathBuf,
    #[arg(short = 'g', long = "genome")]
    genome: PathBuf,
    /// Shift offset (positive = right, negative = left).
    #[arg(short = 's', long = "shift")]
    offset: i64,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        shift(&self.input, &self.genome, self.offset, &mut out)?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Shift BED coordinates by a fixed offset.",
    origin: Some(Origin {
        upstream: "bedtools shift",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["-i INPUT.bed -g genome.txt -s OFFSET"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('i'),
                long: "input",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Input BED.",
                why_default: None,
            },
            FlagSpec {
                short: Some('g'),
                long: "genome",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Genome file.",
                why_default: None,
            },
            FlagSpec {
                short: Some('s'),
                long: "shift",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("i64"),
                required: true,
                default: None,
                description: "Shift offset (+right, -left).",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Shift peaks 500bp downstream",
        command: "rsomics-bed-shift -i peaks.bed -g hg38.genome -s 500",
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
