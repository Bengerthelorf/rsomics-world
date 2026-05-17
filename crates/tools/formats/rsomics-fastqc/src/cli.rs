use std::fs;
use std::io::Write as _;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Context, Result, RsomicsError, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_fastqc::analyze;
use rsomics_fastqc::report::{fastqc_data, summary};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastqc", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// FASTQ file(s). Gzip / bzip2 / xz / zstd auto-detected.
    #[arg(required = true, num_args = 1..)]
    inputs: Vec<PathBuf>,

    /// Output directory; a `<basename>_fastqc/` dir with
    /// `fastqc_data.txt` + `summary.txt` is written per input (`FastQC`
    /// `--extract` layout — `MultiQC` reads `fastqc_data.txt`).
    #[arg(short = 'o', long = "outdir", default_value = ".")]
    outdir: PathBuf,

    /// Write `fastqc_data.txt` for the single input to stdout instead of
    /// to a directory (for piping / `MultiQC` ingestion).
    #[arg(long = "stdout")]
    to_stdout: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(&self) -> Result<()> {
        if self.to_stdout && self.inputs.len() != 1 {
            return Err(RsomicsError::InvalidInput(
                "--stdout requires exactly one input".into(),
            ));
        }
        for input in &self.inputs {
            let mods =
                analyze(input).rs_with_context(|| format!("analyzing {}", input.display()))?;
            let base = input.file_name().map_or_else(
                || input.display().to_string(),
                |n| n.to_string_lossy().into(),
            );
            if self.to_stdout {
                print!("{}", fastqc_data(&mods));
            } else {
                let dir = self.outdir.join(format!("{base}_fastqc"));
                fs::create_dir_all(&dir).map_err(RsomicsError::Io)?;
                fs::File::create(dir.join("fastqc_data.txt"))
                    .map_err(RsomicsError::Io)?
                    .write_all(fastqc_data(&mods).as_bytes())
                    .map_err(RsomicsError::Io)?;
                fs::File::create(dir.join("summary.txt"))
                    .map_err(RsomicsError::Io)?
                    .write_all(summary(&mods, &base).as_bytes())
                    .map_err(RsomicsError::Io)?;
            }
        }
        Ok(())
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Per-file FASTQ quality-control report (independent Rust FastQC reimplementation).",
    origin: Some(Origin {
        upstream: "FastQC",
        upstream_license: "GPL-3.0",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["[OPTIONS] <INPUTS>..."],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('o'),
                long: "outdir",
                aliases: &[],
                value: Some("<DIR>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("."),
                description: "Output dir (writes <base>_fastqc/fastqc_data.txt + summary.txt)",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "stdout",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Emit fastqc_data.txt for one input to stdout",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "json",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Emit AI-friendly JSON envelope on stdout",
                why_default: None,
            },
            FlagSpec {
                short: Some('t'),
                long: "threads",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: None,
                description: "Worker thread count (default: available cores)",
                why_default: None,
            },
            FlagSpec {
                short: Some('h'),
                long: "help",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Show this help (add --plain or --json for alt modes)",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "QC report dir",
            command: "rsomics-fastqc -o qc reads.fastq.gz",
        },
        Example {
            description: "Stream fastqc_data.txt to MultiQC",
            command: "rsomics-fastqc --stdout reads.fq",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    // clap's debug_assert only fires in binary parse, not lib tests — this test is the only way to catch CLI definition errors
    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
