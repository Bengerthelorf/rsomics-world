use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_fastq_sort::{SortKey, sort};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-sort", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(value_name = "INPUT")]
    input: PathBuf,
    #[arg(short = 'l', long = "by-length")]
    by_length: bool,
    #[arg(short = 'L', long = "by-length-desc")]
    by_length_desc: bool,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let key = if self.by_length_desc {
            SortKey::LengthDesc
        } else if self.by_length {
            SortKey::Length
        } else {
            SortKey::Name
        };
        let mut out = std::io::stdout().lock();
        sort(&self.input, key, &mut out)?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Sort FASTQ reads by name or length.",
    origin: Some(Origin {
        upstream: "seqkit sort",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["[-l|-L] <INPUT.fq>"],
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
                description: "Input FASTQ.",
                why_default: None,
            },
            FlagSpec {
                short: Some('l'),
                long: "by-length",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Sort by length.",
                why_default: None,
            },
            FlagSpec {
                short: Some('L'),
                long: "by-length-desc",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Sort longest first.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Sort by name",
        command: "rsomics-fastq-sort reads.fq.gz > sorted.fq",
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
