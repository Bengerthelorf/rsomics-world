use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_gff_stats::stats;
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-gff-stats", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(value_name = "INPUT")]
    input: PathBuf,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let s = stats(&self.input)?;
        println!("Total features:\t{}", s.total);
        println!("\nBy type:");
        for (k, v) in &s.by_type {
            println!("  {k}\t{v}");
        }
        println!("\nBy source:");
        for (k, v) in &s.by_source {
            println!("  {k}\t{v}");
        }
        println!("\nBy chromosome:");
        for (k, v) in &s.by_chrom {
            println!("  {k}\t{v}");
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "GFF/GTF feature statistics by type/source/chromosome.",
    origin: Some(Origin {
        upstream: "AGAT / awk on GFF",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["<INPUT.gff>"],
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
            description: "Input GFF/GTF.",
            why_default: None,
        }],
    }],
    examples: &[Example {
        description: "Feature summary",
        command: "rsomics-gff-stats genes.gff3",
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
