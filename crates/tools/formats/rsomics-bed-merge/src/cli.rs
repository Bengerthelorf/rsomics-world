use std::fs::File;
use std::io::{self, BufWriter};
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use rsomics_intervals::bed;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-merge", disable_help_flag = true)]
pub struct Cli {
    #[arg(short = 'i', long = "input", default_value = "-")]
    input: String,
    #[arg(short = 'o', long = "output", default_value = "-")]
    output: String,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(&self) -> Result<()> {
        let w: BufWriter<Box<dyn io::Write>> = BufWriter::new(if self.output == "-" {
            Box::new(io::stdout().lock())
        } else {
            Box::new(File::create(PathBuf::from(&self.output)).map_err(RsomicsError::Io)?)
        });
        if self.input == "-" {
            bed::merge_sorted(io::stdin().lock(), w)
        } else {
            let f = File::open(PathBuf::from(&self.input)).map_err(RsomicsError::Io)?;
            bed::merge_sorted(f, w)
        }
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Merge overlapping/touching intervals in a pre-sorted BED (bedtools merge equivalent).",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["[OPTIONS] -i <BED>", "[OPTIONS] < input.bed"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('i'),
                long: "input",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("-"),
                description: "Input BED, pre-sorted by chrom then start (bedtools merge contract; pipe rsomics-bed-sort first). Default stdin",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("-"),
                description: "Output BED (default stdout)",
                why_default: None,
            },
            FlagSpec {
                short: Some('h'),
                long: "help",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Show this help (add --plain or --json for alt modes)",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Merge a BED file",
            command: "rsomics-bed-merge -i in.bed -o merged.bed",
        },
        Example {
            description: "Stream from stdin",
            command: "sort -k1,1 -k2,2n in.bed | rsomics-bed-merge > merged.bed",
        },
    ],
    json_result_schema_doc: None,
};
