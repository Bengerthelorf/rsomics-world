use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter};
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use rsomics_intervals::{IntervalSet, bed, complement};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-complement", disable_help_flag = true)]
pub struct Cli {
    #[arg(short = 'i', long = "input")]
    input: PathBuf,
    #[arg(short = 'g', long = "genome")]
    genome: PathBuf,
    #[arg(short = 'o', long = "output", default_value = "-")]
    output: String,
    #[command(flatten)]
    pub common: CommonFlags,
}

fn read_genome(path: &PathBuf) -> Result<Vec<(String, u64)>> {
    let f = File::open(path).map_err(RsomicsError::Io)?;
    let mut out = Vec::new();
    for (lineno, line) in BufReader::new(f).lines().enumerate() {
        let line = line.map_err(RsomicsError::Io)?;
        let trimmed = line.trim_end();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut fields = trimmed.split('\t');
        let chrom = fields.next().ok_or_else(|| {
            RsomicsError::InvalidInput(format!("genome line {}: missing chrom", lineno + 1))
        })?;
        let size_s = fields.next().ok_or_else(|| {
            RsomicsError::InvalidInput(format!("genome line {}: missing size", lineno + 1))
        })?;
        let size: u64 = size_s.parse().map_err(|_| {
            RsomicsError::InvalidInput(format!("genome line {}: bad size {size_s:?}", lineno + 1))
        })?;
        out.push((chrom.to_string(), size));
    }
    Ok(out)
}

impl Cli {
    pub fn execute(&self) -> Result<()> {
        let intervals = bed::read(File::open(&self.input).map_err(RsomicsError::Io)?)?;
        let chrom_sizes = read_genome(&self.genome)?;
        let set: IntervalSet = intervals.into_iter().collect();
        let out = complement(&set, &chrom_sizes);
        let writer: Box<dyn io::Write> = if self.output == "-" {
            Box::new(BufWriter::new(io::stdout().lock()))
        } else {
            Box::new(BufWriter::new(
                File::create(&self.output).map_err(RsomicsError::Io)?,
            ))
        };
        bed::write_bed3(writer, out.iter().cloned())
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Emit regions of each chromosome NOT covered by input BED (bedtools complement equivalent).",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["[OPTIONS] -i <BED> -g <GENOME>"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('i'),
                long: "input",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Input BED",
                why_default: None,
            },
            FlagSpec {
                short: Some('g'),
                long: "genome",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Genome sizes file (TSV: chrom\\tsize)",
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
    examples: &[Example {
        description: "Find gaps in the coverage of a BED",
        command: "rsomics-bed-complement -i covered.bed -g hg38.sizes -o gaps.bed",
    }],
    json_result_schema_doc: None,
};
