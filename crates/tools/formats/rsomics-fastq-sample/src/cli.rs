use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_fastq_sample::sample;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-sample", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Input FASTQ file (gz auto-detected).
    #[arg(value_name = "INPUT")]
    input: PathBuf,

    /// Fraction of records to keep (0.0–1.0).
    #[arg(short = 'p', long = "proportion", default_value_t = 0.1)]
    proportion: f64,

    /// Random seed for reproducibility.
    rng_seed: u64,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        let n = sample(
            &self.input,
            self.proportion,
            self.common.seed_rng(),
            &mut out,
        )?;
        if !self.common.json {
            eprintln!("{n} records sampled");
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Random subsample of FASTQ records.",
    origin: Some(Origin {
        upstream: "seqkit sample",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["[-p FRAC] [--seed N] <INPUT.fq>"],
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
                description: "Input FASTQ (gz auto-detected).",
                why_default: None,
            },
            FlagSpec {
                short: Some('p'),
                long: "proportion",
                aliases: &[],
                value: Some("<F>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("0.1"),
                description: "Fraction of records to keep.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "seed",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("u64"),
                required: false,
                default: Some("42"),
                description: "Random seed.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "10% subsample",
        command: "rsomics-fastq-sample -p 0.1 reads.fq.gz > subset.fq",
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
