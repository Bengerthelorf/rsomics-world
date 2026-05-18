use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bam_sort::{SortOrder, sort_bam};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bam-sort", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Input BAM file.
    #[arg(value_name = "INPUT")]
    input: PathBuf,

    /// Output BAM file.
    #[arg(short = 'o', long = "output")]
    output: PathBuf,

    /// Sort by read name instead of coordinate.
    #[arg(short = 'n', long = "name")]
    by_name: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let order = if self.by_name {
            SortOrder::Name
        } else {
            SortOrder::Coordinate
        };
        sort_bam(&self.input, &self.output, order)
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Sort BAM by coordinate (default) or read name.",
    origin: Some(Origin {
        upstream: "samtools sort",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btp352"),
    }),
    usage_lines: &[
        "<INPUT.bam> -o <OUTPUT.bam>",
        "-n <INPUT.bam> -o <OUTPUT.bam>",
    ],
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
                description: "Input BAM file.",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Output sorted BAM file.",
                why_default: None,
            },
            FlagSpec {
                short: Some('n'),
                long: "name",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Sort by read name (queryname) instead of coordinate.",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Coordinate sort",
            command: "rsomics-bam-sort input.bam -o sorted.bam",
        },
        Example {
            description: "Name sort",
            command: "rsomics-bam-sort -n input.bam -o namesorted.bam",
        },
    ],
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
