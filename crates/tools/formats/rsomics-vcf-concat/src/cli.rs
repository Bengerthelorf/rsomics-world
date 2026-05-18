use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use rsomics_vcf_concat::concat;
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-vcf-concat", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Input VCF files to concatenate.
    #[arg(value_name = "INPUT", required = true)]
    inputs: Vec<PathBuf>,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let paths: Vec<&std::path::Path> = self.inputs.iter().map(PathBuf::as_path).collect();
        let mut out = std::io::stdout().lock();
        let n = concat(&paths, &mut out)?;
        eprintln!("{n} variants concatenated from {} files", self.inputs.len());
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Concatenate VCF files.",
    origin: Some(Origin {
        upstream: "bcftools concat",
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
            description: "VCF files to concatenate.",
            why_default: None,
        }],
    }],
    examples: &[Example {
        description: "Merge per-chromosome VCFs",
        command: "rsomics-vcf-concat chr1.vcf chr2.vcf chr3.vcf > all.vcf",
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
