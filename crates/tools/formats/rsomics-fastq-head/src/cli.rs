use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_fastq_head::head;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-head", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Input FASTQ file (gz auto-detected).
    #[arg(value_name = "INPUT")]
    input: PathBuf,

    /// Number of records to output.
    #[arg(short = 'n', long = "num", default_value_t = 10)]
    num: u64,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        head(&self.input, self.num, &mut out)?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Output the first N FASTQ records.",
    origin: Some(Origin {
        upstream: "seqkit head",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["[-n NUM] <INPUT.fq>"],
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
                short: Some('n'),
                long: "num",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("u64"),
                required: false,
                default: Some("10"),
                description: "Number of records to output.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "First 100 reads",
        command: "rsomics-fastq-head -n 100 reads.fq.gz",
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
