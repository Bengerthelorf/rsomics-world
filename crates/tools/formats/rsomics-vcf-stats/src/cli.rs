use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use rsomics_vcf_stats::stats;
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-vcf-stats", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(value_name = "INPUT")]
    input: PathBuf,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let s = stats(&self.input)?;
        if self.common.json {
            println!("{}", serde_json::to_string_pretty(&s).unwrap_or_default());
        } else {
            println!("Total variants:\t{}", s.total);
            println!("SNPs:\t{}", s.snps);
            println!("Insertions:\t{}", s.insertions);
            println!("Deletions:\t{}", s.deletions);
            println!("MNPs:\t{}", s.mnps);
            println!("Transitions:\t{}", s.transitions);
            println!("Transversions:\t{}", s.transversions);
            println!("Ti/Tv ratio:\t{:.4}", s.ti_tv);
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Basic VCF variant statistics (SNP/indel, Ti/Tv).",
    origin: Some(Origin {
        upstream: "bcftools stats (subset)",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["<INPUT.vcf>"],
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
            description: "Input VCF.",
            why_default: None,
        }],
    }],
    examples: &[Example {
        description: "Variant summary",
        command: "rsomics-vcf-stats variants.vcf",
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
