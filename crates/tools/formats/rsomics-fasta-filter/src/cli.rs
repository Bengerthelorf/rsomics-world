use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_fasta_filter::filter;
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fasta-filter", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(value_name = "INPUT")]
    input: PathBuf,
    /// Minimum sequence length.
    #[arg(short = 'm', long = "min-len", default_value_t = 0)]
    min_len: usize,
    /// Maximum sequence length (0 = no limit).
    #[arg(short = 'M', long = "max-len", default_value_t = 0)]
    max_len: usize,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        filter(&self.input, self.min_len, self.max_len, &mut out)?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Filter FASTA sequences by length.",
    origin: Some(Origin {
        upstream: "seqkit seq -m/-M",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["[-m MIN] [-M MAX] <INPUT.fa>"],
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
                short: Some('m'),
                long: "min-len",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("0"),
                description: "Min length.",
                why_default: None,
            },
            FlagSpec {
                short: Some('M'),
                long: "max-len",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("0"),
                description: "Max length (0=no limit).",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Keep contigs >= 1kb",
        command: "rsomics-fasta-filter -m 1000 assembly.fa > long.fa",
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
