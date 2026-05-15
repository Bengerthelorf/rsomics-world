use std::fs::File;
use std::io::{self, BufWriter};
use std::path::PathBuf;
use std::process::{self, ExitCode};

use clap::Parser;
use rsomics_common::{CommonFlags, ExitCode as RsomicsExit, Result, RsomicsError, ToolMeta, run};
use rsomics_help::{
    Example, FlagSpec, HelpSpec, Origin, Section, intercept_help, render as render_help,
};
use rsomics_intervals::bed;

const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

const HELP: HelpSpec = HelpSpec {
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

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-merge", disable_help_flag = true)]
struct Cli {
    #[arg(short = 'i', long = "input", default_value = "-")]
    input: String,
    #[arg(short = 'o', long = "output", default_value = "-")]
    output: String,
    #[command(flatten)]
    common: CommonFlags,
}

fn pipeline(cli: &Cli) -> Result<()> {
    let w: BufWriter<Box<dyn io::Write>> = BufWriter::new(if cli.output == "-" {
        Box::new(io::stdout().lock())
    } else {
        Box::new(File::create(PathBuf::from(&cli.output)).map_err(RsomicsError::Io)?)
    });
    if cli.input == "-" {
        bed::merge_sorted(io::stdin().lock(), w)
    } else {
        let f = File::open(PathBuf::from(&cli.input)).map_err(RsomicsError::Io)?;
        bed::merge_sorted(f, w)
    }
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
