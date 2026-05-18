use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bam_index::index_bam;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bam-index", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Input coordinate-sorted BAM file.
    #[arg(value_name = "INPUT")]
    input: PathBuf,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        index_bam(&self.input)
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Create BAI index for a coordinate-sorted BAM.",
    origin: Some(Origin {
        upstream: "samtools index",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btp352"),
    }),
    usage_lines: &["<INPUT.bam>"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[FlagSpec {
            short: None,
            long: "INPUT",
            aliases: &[],
            value: Some("<path>"),
            type_hint: Some("Path"),
            required: true,
            default: None,
            description: "Coordinate-sorted BAM file to index.",
            why_default: None,
        }],
    }],
    examples: &[Example {
        description: "Index a sorted BAM",
        command: "rsomics-bam-index sorted.bam",
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
