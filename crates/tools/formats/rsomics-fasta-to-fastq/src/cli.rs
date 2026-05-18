use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_fasta_to_fastq::convert;
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fasta-to-fastq", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(value_name = "INPUT")]
    input: PathBuf,
    /// Phred quality score to assign (default 40 = 'I').
    #[arg(long = "qual", default_value_t = 40)]
    qual: u8,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        convert(&self.input, self.qual + 33, &mut out)?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Convert FASTA to FASTQ with constant quality.",
    origin: Some(Origin {
        upstream: "seqkit fq2fa (reverse)",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["[--qual Q] <INPUT.fa>"],
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
                long: "qual",
                aliases: &[],
                value: Some("<Q>"),
                type_hint: Some("u8"),
                required: false,
                default: Some("40"),
                description: "Phred quality to assign.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "FASTA to FASTQ with Q40",
        command: "rsomics-fasta-to-fastq ref.fa > ref.fq",
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
