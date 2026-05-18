use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_fasta_split::split_by_count;
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fasta-split", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(value_name = "INPUT")]
    input: PathBuf,
    /// Sequences per output file.
    #[arg(long = "seqs-per-file", default_value_t = 1000)]
    seqs_per_file: u64,
    /// Output file prefix.
    #[arg(long = "prefix", default_value = "split_")]
    prefix: String,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let n = split_by_count(&self.input, self.seqs_per_file, &self.prefix)?;
        eprintln!("{n} files written");
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Split FASTA into multiple files by sequence count.",
    origin: Some(Origin {
        upstream: "seqkit split2 / pyfasta split",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["[--seqs-per-file N] <INPUT.fa>"],
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
                short: None,
                long: "seqs-per-file",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("u64"),
                required: false,
                default: Some("1000"),
                description: "Sequences per file.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "prefix",
                aliases: &[],
                value: Some("<str>"),
                type_hint: Some("String"),
                required: false,
                default: Some("split_"),
                description: "Output prefix.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Split genome into 100-seq chunks",
        command: "rsomics-fasta-split --seqs-per-file 100 genome.fa",
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
