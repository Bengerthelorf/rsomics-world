use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bam_idxstats::idxstats;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bam-idxstats", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Input BAM file (must have a .bai index).
    #[arg(value_name = "INPUT")]
    input: PathBuf,

    /// Output format: `default` (samtools-compatible TSV) or `json`.
    #[arg(short = 'O', long = "output-fmt", default_value = "default")]
    output_fmt: OutputFmt,

    #[command(flatten)]
    pub common: CommonFlags,
}

#[derive(Debug, Clone, clap::ValueEnum)]
enum OutputFmt {
    Default,
    Json,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let stats = idxstats(&self.input)?;

        match self.output_fmt {
            OutputFmt::Default => print!("{stats}"),
            OutputFmt::Json => {
                let json = serde_json::to_string_pretty(&stats)
                    .map_err(|e| RsomicsError::InvalidInput(format!("json: {e}")))?;
                println!("{json}");
            }
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Per-reference mapped/unmapped read counts from a BAM index.",
    origin: Some(Origin {
        upstream: "samtools idxstats",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btp352"),
    }),
    usage_lines: &["<INPUT.bam>"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[FlagSpec {
            short: None,
            long: "INPUT",
            aliases: &[],
            value: Some("<path>"),
            type_hint: Some("Path"),
            required: true,
            default: None,
            description: "BAM file with .bai index.",
            why_default: None,
        }],
    }],
    examples: &[Example {
        description: "Print per-reference read counts",
        command: "rsomics-bam-idxstats aligned.bam",
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
