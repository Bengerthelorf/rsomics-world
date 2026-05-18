use clap::Parser;
use rsomics_bed_random::random_bed;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-random", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Genome file (chrom<TAB>size).
    #[arg(short = 'g', long = "genome")]
    genome: PathBuf,
    /// Number of intervals to generate.
    #[arg(short = 'n', long = "num", default_value_t = 1000)]
    num: u64,
    /// Length of each interval.
    #[arg(short = 'l', long = "length", default_value_t = 1000)]
    length: u64,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        random_bed(
            &self.genome,
            self.num,
            self.length,
            self.common.seed_rng(),
            &mut out,
        )
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Generate random BED intervals.",
    origin: Some(Origin {
        upstream: "bedtools random",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["-g genome.txt [-n NUM] [-l LEN]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
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
                short: Some('n'),
                long: "num",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("u64"),
                required: false,
                default: Some("1000"),
                description: "Number of intervals.",
                why_default: None,
            },
            FlagSpec {
                short: Some('l'),
                long: "length",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("u64"),
                required: false,
                default: Some("1000"),
                description: "Interval length.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "10k random 500bp regions",
        command: "rsomics-bed-random -g hg38.genome -n 10000 -l 500 --seed 42",
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
