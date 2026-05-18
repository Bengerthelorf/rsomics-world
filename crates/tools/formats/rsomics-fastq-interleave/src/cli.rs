use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_fastq_interleave::interleave;
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-interleave", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// R1 FASTQ file.
    #[arg(long = "in1", short = 'i')]
    in1: PathBuf,
    /// R2 FASTQ file.
    #[arg(long = "in2", short = 'I')]
    in2: PathBuf,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        interleave(&self.in1, &self.in2, &mut out)?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Interleave paired FASTQ R1/R2 into a single stream.",
    origin: Some(Origin {
        upstream: "seqkit pair / BBMap repair.sh",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["--in1 R1.fq --in2 R2.fq > interleaved.fq"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('i'),
                long: "in1",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "R1 FASTQ.",
                why_default: None,
            },
            FlagSpec {
                short: Some('I'),
                long: "in2",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "R2 FASTQ.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Interleave for BWA",
        command: "rsomics-fastq-interleave -i R1.fq.gz -I R2.fq.gz | bwa mem -p ref.fa -",
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
