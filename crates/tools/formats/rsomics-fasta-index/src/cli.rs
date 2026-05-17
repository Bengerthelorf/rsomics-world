use std::path::PathBuf;

use clap::{Parser, Subcommand};
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_fasta_index::{build_index, fetch_region, write_index};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fasta-index", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[command(subcommand)]
    cmd: Cmd,

    #[command(flatten)]
    pub common: CommonFlags,
}

#[derive(Subcommand, Debug)]
enum Cmd {
    /// Create a .fai index for a FASTA file.
    Index {
        /// Input FASTA file.
        #[arg(value_name = "FASTA")]
        input: PathBuf,
    },
    /// Fetch a region from an indexed FASTA.
    Fetch {
        /// Input FASTA file (must have .fai index).
        #[arg(value_name = "FASTA")]
        input: PathBuf,
        /// Region(s) to fetch: `chr1`, `chr1:100-200` (1-based, inclusive).
        #[arg(value_name = "REGION")]
        regions: Vec<String>,
    },
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        match self.cmd {
            Cmd::Index { input } => {
                let index = build_index(&input)?;
                let fai_path = input.with_extension(format!(
                    "{}.fai",
                    input.extension().unwrap_or_default().to_string_lossy()
                ));
                write_index(&index, &fai_path)?;
                eprintln!("wrote {}", fai_path.display());
                Ok(())
            }
            Cmd::Fetch { input, regions } => {
                if regions.is_empty() {
                    return Err(RsomicsError::InvalidInput("no region specified".into()));
                }
                let fai_path = input.with_extension(format!(
                    "{}.fai",
                    input.extension().unwrap_or_default().to_string_lossy()
                ));
                for region in &regions {
                    let seq = fetch_region(&input, &fai_path, region)?;
                    println!(">{region}");
                    for chunk in seq.chunks(80) {
                        println!("{}", std::str::from_utf8(chunk).unwrap_or("?"));
                    }
                }
                Ok(())
            }
        }
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "FASTA index (.fai) creation and random-access fetch.",
    origin: Some(Origin {
        upstream: "samtools faidx",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btp352"),
    }),
    usage_lines: &["index <FASTA>", "fetch <FASTA> <REGION>..."],
    sections: &[Section {
        title: "SUBCOMMANDS",
        flags: &[
            FlagSpec {
                short: None,
                long: "index",
                aliases: &[],
                value: Some("<FASTA>"),
                type_hint: Some("Path"),
                required: false,
                default: None,
                description: "Create .fai index for a FASTA file.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "fetch",
                aliases: &[],
                value: Some("<FASTA> <REGION>..."),
                type_hint: Some("Path + regions"),
                required: false,
                default: None,
                description: "Fetch subsequences by region from indexed FASTA.",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Index a FASTA",
            command: "rsomics-fasta-index index ref.fa",
        },
        Example {
            description: "Fetch chr1:1000-2000",
            command: "rsomics-fasta-index fetch ref.fa chr1:1000-2000",
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
