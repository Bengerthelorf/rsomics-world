use std::io::BufReader;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bam_flagstat::{count_bam, count_sam};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bam-flagstat", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Input BAM/SAM file. Use `-` for stdin (SAM only).
    #[arg(value_name = "INPUT")]
    input: PathBuf,

    /// Output format: `default` (samtools-compatible text) or `json`.
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
        let counts = if self.input.as_os_str() == "-" {
            let stdin = std::io::stdin().lock();
            count_sam(BufReader::new(stdin))?
        } else {
            let ext = self.input.extension().and_then(|e| e.to_str()).unwrap_or("");
            match ext {
                "sam" => {
                    let file = std::fs::File::open(&self.input)
                        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", self.input.display())))?;
                    count_sam(BufReader::new(file))?
                }
                _ => count_bam(&self.input)?,
            }
        };

        match self.output_fmt {
            OutputFmt::Default => print!("{counts}"),
            OutputFmt::Json => {
                let json = serde_json::to_string_pretty(&counts)
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
    tagline: "SAM/BAM flag statistics (per-flag read counts, QC-pass/fail split).",
    origin: Some(Origin {
        upstream: "samtools flagstat",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btp352"),
    }),
    usage_lines: &[
        "<INPUT.bam>",
        "-O json <INPUT.bam>",
    ],
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
                description: "BAM, SAM, or CRAM file (positional). `-` reads SAM from stdin.",
                why_default: None,
            },
            FlagSpec {
                short: Some('O'),
                long: "output-fmt",
                aliases: &[],
                value: Some("<fmt>"),
                type_hint: Some("String"),
                required: false,
                default: Some("default"),
                description: "Output format: `default` (samtools text) or `json`.",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Print flag statistics",
            command: "rsomics-bam-flagstat aligned.bam",
        },
        Example {
            description: "JSON output",
            command: "rsomics-bam-flagstat -O json aligned.bam",
        },
    ],
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
