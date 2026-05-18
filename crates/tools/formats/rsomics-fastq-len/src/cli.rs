use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_fastq_len::lengths;
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-len", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(value_name = "INPUT")]
    input: PathBuf,
    #[arg(long = "tab")]
    tab: bool,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        lengths(&self.input, self.tab, &mut out)?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Output per-read lengths from FASTQ.",
    origin: Some(Origin {
        upstream: "seqkit fx2tab -nl",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["[--tab] <INPUT.fq>"],
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
                long: "tab",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Output name<TAB>length.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Read length distribution",
        command: "rsomics-fastq-len reads.fq.gz | sort -n | uniq -c",
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
