use clap::Parser;
use rsomics_bed_makewindows::makewindows;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-makewindows", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Genome file (chrom<TAB>size).
    #[arg(short = 'g', long = "genome")]
    genome: PathBuf,
    /// Window size in bases.
    #[arg(short = 'w', long = "window", default_value_t = 1000)]
    window: u64,
    /// Step size (default = window size, i.e. non-overlapping).
    #[arg(short = 's', long = "step")]
    step: Option<u64>,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let step = self.step.unwrap_or(self.window);
        let mut out = std::io::stdout().lock();
        makewindows(&self.genome, self.window, step, &mut out)?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Divide genome into fixed-size tiling windows.",
    origin: Some(Origin {
        upstream: "bedtools makewindows",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["-g genome.txt -w SIZE [-s STEP]"],
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
                short: Some('w'),
                long: "window",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("u64"),
                required: false,
                default: Some("1000"),
                description: "Window size.",
                why_default: None,
            },
            FlagSpec {
                short: Some('s'),
                long: "step",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("u64"),
                required: false,
                default: None,
                description: "Step size (default=window).",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "1kb non-overlapping tiles",
        command: "rsomics-bed-makewindows -g hg38.genome -w 1000",
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
