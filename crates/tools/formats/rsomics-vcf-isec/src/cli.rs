use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use rsomics_vcf_isec::isec;
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-vcf-isec", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// First VCF (variants to filter).
    #[arg(short = 'a')]
    file_a: PathBuf,
    /// Second VCF (reference set).
    #[arg(short = 'b')]
    file_b: PathBuf,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let mut out = std::io::stdout().lock();
        let n = isec(&self.file_a, &self.file_b, &mut out)?;
        eprintln!("{n} shared variants");
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Find shared variants between two VCFs.",
    origin: Some(Origin {
        upstream: "bcftools isec",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["-a A.vcf -b B.vcf"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('a'),
                long: "a",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "VCF A.",
                why_default: None,
            },
            FlagSpec {
                short: Some('b'),
                long: "b",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "VCF B (reference set).",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Shared variants",
        command: "rsomics-vcf-isec -a sample.vcf -b dbsnp.vcf > shared.vcf",
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
