use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_fastq_rename::rename;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-rename", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Input FASTQ file (gz auto-detected).
    #[arg(value_name = "INPUT")]
    input: PathBuf,

    /// Prefix for new read names (sequential: prefix0, prefix1, ...).
    #[arg(long = "prefix", default_value = "read_")]
    prefix: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        rename(&self.input, &self.prefix, &mut out)?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Rename FASTQ reads with sequential IDs.",
    origin: Some(Origin {
        upstream: "seqkit rename",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["[--prefix STR] <INPUT.fq>"],
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
                short: None,
                long: "prefix",
                aliases: &[],
                value: Some("<str>"),
                type_hint: Some("String"),
                required: false,
                default: Some("read_"),
                description: "Prefix for sequential IDs.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Rename with prefix",
        command: "rsomics-fastq-rename --prefix sample1_ reads.fq > renamed.fq",
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
