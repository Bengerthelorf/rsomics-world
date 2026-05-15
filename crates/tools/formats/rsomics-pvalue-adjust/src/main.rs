use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::process::{self, ExitCode};

use clap::Parser;
use rsomics_common::{CommonFlags, ExitCode as RsomicsExit, Result, RsomicsError, ToolMeta, run};
use rsomics_help::{
    Example, FlagSpec, HelpSpec, Origin, Section, intercept_help, render as render_help,
};
use rsomics_stats::{bh_adjust, bonferroni_adjust};

const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Apply BH FDR or Bonferroni correction to a column of p-values (R's `p.adjust` equivalent).",
    origin: Some(Origin {
        upstream: "R stats::p.adjust",
        upstream_license: "GPL",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1111/j.2517-6161.1995.tb02031.x"),
    }),
    usage_lines: &[
        "[OPTIONS] -i <FILE>",
        "cat pvals.txt | rsomics-pvalue-adjust",
    ],
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
                description: "Input file (one p-value per line; `-` = stdin)",
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
                description: "Output file (tab-separated: original\\tadjusted; `-` = stdout)",
                why_default: None,
            },
            FlagSpec {
                short: Some('m'),
                long: "method",
                aliases: &[],
                value: Some("<bh|bonferroni>"),
                type_hint: Some("enum"),
                required: false,
                default: Some("bh"),
                description: "Correction method (bh = Benjamini-Hochberg FDR, bonferroni = family-wise)",
                why_default: None,
            },
            FlagSpec {
                short: Some('c'),
                long: "column",
                aliases: &[],
                value: Some("<n>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("1"),
                description: "1-based column index when input has multiple tab-separated columns",
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
            description: "BH-adjust a single column from stdin",
            command: "cut -f5 results.tsv | rsomics-pvalue-adjust > adjusted.tsv",
        },
        Example {
            description: "Bonferroni on column 7 of a TSV",
            command: "rsomics-pvalue-adjust -i de.tsv -c 7 -m bonferroni -o adjusted.tsv",
        },
    ],
    json_result_schema_doc: None,
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-pvalue-adjust", disable_help_flag = true)]
struct Cli {
    #[arg(short = 'i', long = "input", default_value = "-")]
    input: String,
    #[arg(short = 'o', long = "output", default_value = "-")]
    output: String,
    #[arg(short = 'm', long = "method", default_value = "bh")]
    method: String,
    #[arg(short = 'c', long = "column", default_value_t = 1)]
    column: usize,
    #[command(flatten)]
    common: CommonFlags,
}

fn pipeline(cli: &Cli) -> Result<()> {
    let lines = read_lines(&cli.input)?;
    let (raw_lines, pvals) = parse_column(&lines, cli.column)?;
    let adjusted = match cli.method.as_str() {
        "bh" => bh_adjust(&pvals),
        "bonferroni" => bonferroni_adjust(&pvals),
        other => {
            return Err(RsomicsError::ConfigError(format!(
                "unknown --method {other:?} (expected `bh` or `bonferroni`)"
            )));
        }
    }
    .map_err(|e| RsomicsError::InvalidInput(format!("p-value adjust: {e}")))?;
    write_out(&cli.output, &raw_lines, &pvals, &adjusted)
}

fn read_lines(path: &str) -> Result<Vec<String>> {
    let reader: Box<dyn BufRead> = if path == "-" {
        Box::new(BufReader::new(io::stdin().lock()))
    } else {
        Box::new(BufReader::new(File::open(path).map_err(RsomicsError::Io)?))
    };
    reader
        .lines()
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(RsomicsError::Io)
}

fn parse_column(lines: &[String], col: usize) -> Result<(Vec<String>, Vec<f64>)> {
    if col == 0 {
        return Err(RsomicsError::ConfigError(
            "--column is 1-based; 0 is not valid".into(),
        ));
    }
    let mut kept = Vec::new();
    let mut pvals = Vec::new();
    for (lineno, line) in lines.iter().enumerate() {
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        }
        let field = trimmed.split('\t').nth(col - 1).ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "line {} has fewer than {col} tab-separated columns",
                lineno + 1
            ))
        })?;
        let p: f64 = field.parse().map_err(|_| {
            RsomicsError::InvalidInput(format!("line {}: bad p-value {field:?}", lineno + 1))
        })?;
        pvals.push(p);
        kept.push(trimmed.to_string());
    }
    Ok((kept, pvals))
}

fn write_out(path: &str, lines: &[String], pvals: &[f64], adjusted: &[f64]) -> Result<()> {
    let mut w: Box<dyn Write> = if path == "-" {
        Box::new(BufWriter::new(io::stdout().lock()))
    } else {
        Box::new(BufWriter::new(
            File::create(path).map_err(RsomicsError::Io)?,
        ))
    };
    for ((line, p), adj) in lines.iter().zip(pvals.iter()).zip(adjusted.iter()) {
        if line.contains('\t') {
            writeln!(w, "{line}\t{adj:.6}").map_err(RsomicsError::Io)?;
        } else {
            writeln!(w, "{p}\t{adj:.6}").map_err(RsomicsError::Io)?;
        }
    }
    Ok(())
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
