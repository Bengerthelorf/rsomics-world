use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use rsomics_vcf_merge::merge;
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-vcf-merge", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(value_name = "INPUT", required = true)]
    inputs: Vec<PathBuf>,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let paths: Vec<&std::path::Path> = self.inputs.iter().map(PathBuf::as_path).collect();
        let mut out = std::io::stdout().lock();
        let n = merge(&paths, &mut out)?;
        eprintln!("{n} variants merged from {} files", self.inputs.len());
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Merge multi-sample VCFs by position.",
    origin: Some(Origin {
        upstream: "bcftools merge",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["<A.vcf> <B.vcf> [C.vcf ...]"],
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
            description: "VCF files to merge.",
            why_default: None,
        }],
    }],
    examples: &[Example {
        description: "Merge per-sample VCFs",
        command: "rsomics-vcf-merge sample1.vcf sample2.vcf > merged.vcf",
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
