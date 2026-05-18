use clap::Parser;
use rsomics_bed_slop::slop;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-slop", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Input BED file.
    #[arg(short = 'i', long = "input")]
    input: PathBuf,
    /// Genome file (chrom<TAB>size).
    #[arg(short = 'g', long = "genome")]
    genome: PathBuf,
    /// Bases to add on the left (5') side.
    #[arg(short = 'l', long = "left", default_value_t = 0)]
    left: u64,
    /// Bases to add on the right (3') side.
    #[arg(short = 'r', long = "right", default_value_t = 0)]
    right: u64,
    /// Bases to add on both sides (overrides -l/-r).
    #[arg(short = 'b', long = "both", default_value_t = 0)]
    both: u64,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let (l, r) = if self.both > 0 {
            (self.both, self.both)
        } else {
            (self.left, self.right)
        };
        let mut out = std::io::stdout().lock();
        slop(&self.input, &self.genome, l, r, &mut out)?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Extend BED intervals by a fixed amount.",
    origin: Some(Origin {
        upstream: "bedtools slop",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["-i INPUT.bed -g genome.txt -b 100"],
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
                description: "Genome file (chrom<TAB>size).",
                why_default: None,
            },
            FlagSpec {
                short: Some('b'),
                long: "both",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("u64"),
                required: false,
                default: Some("0"),
                description: "Extend both sides.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Extend 500bp each side",
        command: "rsomics-bed-slop -i peaks.bed -g hg38.genome -b 500",
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
