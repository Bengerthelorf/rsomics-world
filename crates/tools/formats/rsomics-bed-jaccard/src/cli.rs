use clap::Parser;
use rsomics_bed_jaccard::jaccard;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-jaccard", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(short = 'a')]
    file_a: PathBuf,
    #[arg(short = 'b')]
    file_b: PathBuf,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let r = jaccard(&self.file_a, &self.file_b)?;
        println!("intersection\tunion\tjaccard\tn_intersections");
        println!(
            "{}\t{}\t{:.6}\t{}",
            r.intersection, r.union, r.jaccard, r.n_intersections
        );
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Compute Jaccard similarity between two BED files.",
    origin: Some(Origin {
        upstream: "bedtools jaccard",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["-a A.bed -b B.bed"],
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
                description: "BED file A.",
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
                description: "BED file B.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Jaccard between peaks and enhancers",
        command: "rsomics-bed-jaccard -a peaks.bed -b enhancers.bed",
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
