use std::io::BufWriter;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_vcf_filter::{FilterConfig, filter_vcf};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-vcf-filter", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Input VCF/BCF file.
    #[arg(value_name = "INPUT")]
    input: PathBuf,

    /// Output file (default stdout).
    #[arg(short = 'o', long = "output", default_value = "-")]
    output: String,

    /// Keep only PASS variants (FILTER == PASS or .).
    #[arg(short = 'f', long = "pass-only")]
    pass_only: bool,

    /// Minimum QUAL score.
    #[arg(short = 'q', long = "min-qual")]
    min_qual: Option<f32>,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let cfg = FilterConfig {
            min_qual: self.min_qual,
            pass_only: self.pass_only,
            ..FilterConfig::default()
        };

        let mut out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(BufWriter::new(std::io::stdout().lock()))
        } else {
            Box::new(BufWriter::new(
                std::fs::File::create(&self.output).map_err(RsomicsError::Io)?,
            ))
        };

        let stats = filter_vcf(&self.input, &mut out, &cfg)?;

        if !self.common.json {
            eprintln!(
                "{}/{} variants passed filter",
                stats.passed, stats.total
            );
        }

        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "VCF record filtering by quality, FILTER status, and expressions.",
    origin: Some(Origin {
        upstream: "bcftools view/filter",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/gigascience/giab008"),
    }),
    usage_lines: &["[OPTIONS] <INPUT.vcf>"],
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
                description: "Input VCF/BCF file.",
                why_default: None,
            },
            FlagSpec {
                short: Some('f'),
                long: "pass-only",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Keep only PASS variants.",
                why_default: None,
            },
            FlagSpec {
                short: Some('q'),
                long: "min-qual",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("f32"),
                required: false,
                default: None,
                description: "Minimum QUAL score.",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Keep only PASS variants with QUAL >= 30",
            command: "rsomics-vcf-filter --pass-only --min-qual 30 input.vcf > filtered.vcf",
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
