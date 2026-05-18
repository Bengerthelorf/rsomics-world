use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_fasta_sample::sample;
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fasta-sample", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(value_name = "INPUT")]
    input: PathBuf,
    #[arg(short = 'p', long = "proportion", default_value_t = 0.1)]
    proportion: f64,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        sample(
            &self.input,
            self.proportion,
            self.common.seed_rng(),
            &mut out,
        )?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Random subsample of FASTA sequences.",
    origin: Some(Origin {
        upstream: "seqkit sample (FASTA)",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["[-p FRAC] <INPUT.fa>"],
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
                description: "Input FASTA.",
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
                description: "Fraction to keep.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "10% subsample",
        command: "rsomics-fasta-sample -p 0.1 genome.fa --seed 42 > subset.fa",
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
