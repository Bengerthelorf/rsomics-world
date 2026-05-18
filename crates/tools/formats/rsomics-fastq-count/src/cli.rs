use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_fastq_count::count;
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-count", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Input FASTQ file(s).
    #[arg(value_name = "INPUT", required = true)]
    inputs: Vec<PathBuf>,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut total: u64 = 0;
        for input in &self.inputs {
            let n = count(input)?;
            if self.inputs.len() > 1 {
                println!("{n}\t{}", input.display());
            }
            total += n;
        }
        if self.inputs.len() == 1 {
            println!("{total}");
        } else {
            println!("{total}\ttotal");
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Count FASTQ records (gz-transparent).",
    origin: Some(Origin {
        upstream: "seqkit stats -T (record count only)",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["<INPUT.fq> [INPUT2.fq ...]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[FlagSpec {
            short: None,
            long: "INPUT",
            aliases: &[],
            value: Some("<path>..."),
            type_hint: Some("Path"),
            required: true,
            default: None,
            description: "Input FASTQ file(s).",
            why_default: None,
        }],
    }],
    examples: &[Example {
        description: "Count reads in a file",
        command: "rsomics-fastq-count reads.fq.gz",
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
