use std::fs::File;
use std::io::{self, BufWriter};
use std::path::PathBuf;
use std::process::{self, ExitCode};

use clap::Parser;
use rsomics_common::{CommonFlags, ExitCode as RsomicsExit, Result, ToolMeta, run};
use rsomics_help::{
    Example, FlagSpec, HelpSpec, Origin, Section, intercept_help, render as render_help,
};
use rsomics_intervals::{IntervalSet, bed, intersect};

const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Intersect two BED files, emitting clipped overlap regions (bedtools intersect equivalent).",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["[OPTIONS] -a <BED> -b <BED>"],
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
                description: "Left input BED",
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
                description: "Right input BED",
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
        description: "Intersect peaks and gene bodies",
        command: "rsomics-bed-intersect -a peaks.bed -b genes.bed -o overlaps.bed",
    }],
    json_result_schema_doc: None,
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-intersect", disable_help_flag = true)]
struct Cli {
    #[arg(short = 'a', long = "a")]
    a: PathBuf,
    #[arg(short = 'b', long = "b")]
    b: PathBuf,
    #[arg(short = 'o', long = "output", default_value = "-")]
    output: String,
    #[command(flatten)]
    common: CommonFlags,
}

fn pipeline(cli: &Cli) -> Result<()> {
    let a_ivs = bed::read(File::open(&cli.a).map_err(rsomics_common::RsomicsError::Io)?)?;
    let b_ivs = bed::read(File::open(&cli.b).map_err(rsomics_common::RsomicsError::Io)?)?;
    let a_set: IntervalSet = a_ivs.into_iter().collect();
    let b_set: IntervalSet = b_ivs.into_iter().collect();
    let out = intersect(&a_set, &b_set);
    let writer: Box<dyn io::Write> = if cli.output == "-" {
        Box::new(BufWriter::new(io::stdout().lock()))
    } else {
        Box::new(BufWriter::new(
            File::create(&cli.output).map_err(rsomics_common::RsomicsError::Io)?,
        ))
    };
    bed::write_bed3(writer, out.iter().cloned())
}

fn main() -> ExitCode {
    let raw_args: Vec<String> = std::env::args().collect();
    if let Some(mode) = intercept_help(&raw_args) {
        render_help(&HELP, mode);
        return process::ExitCode::from(RsomicsExit::Ok);
    }
    let cli = Cli::parse();
    let common = cli.common.clone();
    run(&common, META, || pipeline(&cli))
}
