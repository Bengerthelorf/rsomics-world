use clap::Parser;
use rsomics_bed_stats::stats;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use std::path::PathBuf;
pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};
#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-stats", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    #[arg(value_name = "INPUT")]
    input: PathBuf,
    #[command(flatten)]
    pub common: CommonFlags,
}
impl Cli {
    pub fn execute(self) -> Result<()> {
        let s = stats(&self.input)?;
        println!("Intervals:\t{}", s.count);
        println!("Total bases:\t{}", s.total_bases);
        println!("Min length:\t{}", s.min_len);
        println!("Max length:\t{}", s.max_len);
        println!("Mean length:\t{:.1}", s.mean_len);
        println!("Median length:\t{:.1}", s.median_len);
        Ok(())
    }
}
pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "BED interval statistics.",
    origin: Some(Origin {
        upstream: "bedtools / awk stats",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["<INPUT.bed>"],
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
            description: "Input BED.",
            why_default: None,
        }],
    }],
    examples: &[Example {
        description: "Interval summary",
        command: "rsomics-bed-stats peaks.bed",
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
