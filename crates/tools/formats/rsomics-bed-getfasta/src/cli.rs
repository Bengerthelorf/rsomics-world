use clap::Parser;
use rsomics_bed_getfasta::getfasta;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-getfasta", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// BED file with regions to extract.
    #[arg(long = "bed")]
    bed: PathBuf,
    /// Indexed FASTA reference (must have .fai).
    #[arg(long = "fi")]
    fasta: PathBuf,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        getfasta(&self.bed, &self.fasta, &mut out)?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Extract FASTA sequences for BED regions from an indexed reference.",
    origin: Some(Origin {
        upstream: "bedtools getfasta",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["--bed regions.bed --fi ref.fa"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "bed",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "BED file with regions.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "fi",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Indexed FASTA (needs .fai).",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Extract peak sequences",
        command: "rsomics-bed-getfasta --bed peaks.bed --fi hg38.fa > peaks.fa",
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
