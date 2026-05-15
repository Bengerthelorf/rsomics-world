use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use rsomics_common::{CommonFlags, Context, Result, RsomicsError, ToolMeta, run};
use rsomics_help::{
    Example, FlagSpec, HelpSpec, Origin, Section, intercept_help, render as render_help,
};

use rsomics_fasta_stats::{Config, FastaStats, compute_stats, render_pretty, render_tabular};

const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

// Declarative help spec. The binary owns the data only — rich / plain /
// json layouts all dispatch from this single struct in `rsomics-help`.
const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Per-file statistics for FASTA inputs (Rust port of seqkit stats — FASTA subset).",
    origin: Some(Origin {
        upstream: "seqkit",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1371/journal.pone.0163962"),
    }),
    usage_lines: &["[OPTIONS] <INPUTS>..."],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('a'),
                long: "all",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Emit extended stats (Q1/Q2/Q3, N50, GC%, sum_gap, sum_n)",
                why_default: Some(
                    "seqkit's `--all` opts in to the extended columns; off by default",
                ),
            },
            FlagSpec {
                short: Some('T'),
                long: "tabular",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Tab-separated machine-readable output",
                why_default: None,
            },
            FlagSpec {
                short: Some('G'),
                long: "gap-letters",
                aliases: &[],
                value: Some("<CHARS>"),
                type_hint: Some("String"),
                required: false,
                default: Some("- ."),
                description: "Characters counted as gap when --all is set",
                why_default: Some("matches seqkit's default gap letter set"),
            },
            FlagSpec {
                short: None,
                long: "json",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Emit AI-friendly JSON envelope on stdout",
                why_default: None,
            },
            FlagSpec {
                short: Some('t'),
                long: "threads",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: None,
                description: "Worker thread count (default: available cores)",
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
            description: "Default stats",
            command: "rsomics-fasta-stats genome.fa",
        },
        Example {
            description: "Extended, tabular, on gzip input",
            command: "rsomics-fasta-stats --tabular --all genome.fa.gz",
        },
        Example {
            description: "JSON envelope through jq",
            command: "rsomics-fasta-stats --json scaffolds.fa | jq .result",
        },
    ],
    json_result_schema_doc: Some("https://docs.rs/rsomics-fasta-stats/0.4/#json-output-schema"),
};

/// Compute per-file statistics for FASTA inputs.
///
/// Output is the FASTA subset of `seqkit stats`. For FASTQ-only quality
/// columns, use the sibling FASTQ-stats crate.
#[derive(Parser, Debug)]
#[command(name = "rsomics-fasta-stats", version, about, long_about = None, disable_help_flag = true)]
struct Cli {
    /// FASTA file(s). Gzip / bzip2 / xz / zstd inputs are auto-detected
    /// by extension or magic bytes. Stdin (`-`) is not yet supported.
    #[arg(required = true, num_args = 1..)]
    inputs: Vec<PathBuf>,

    /// Emit extended statistics: `Q1` / `Q2` / `Q3`, `sum_gap`, `N50`,
    /// `N50_num`, `GC(%)`, `sum_n`.
    #[arg(short = 'a', long = "all")]
    all: bool,

    /// Tab-separated machine-readable output. One header line, one row per
    /// input. Disables thousand-separator commas.
    #[arg(short = 'T', long = "tabular")]
    tabular: bool,

    /// Characters counted as gap when `--all` is set. Default matches
    /// seqkit: hyphen, space, period.
    #[arg(short = 'G', long = "gap-letters", default_value = "- .")]
    gap_letters: String,

    #[command(flatten)]
    common: CommonFlags,
}

fn pipeline(args: &Cli) -> Result<Vec<FastaStats>> {
    let cfg = Config {
        extended: args.all,
        gap_letters: args.gap_letters.as_bytes().to_vec(),
    };
    let mut results = Vec::with_capacity(args.inputs.len());
    for input in &args.inputs {
        if input.as_os_str() == "-" {
            return Err(RsomicsError::InvalidInput(
                "stdin (`-`) input not yet supported; pass a file path".into(),
            ));
        }
        let stats = compute_stats(input, &cfg)
            .rs_with_context(|| format!("computing stats for {}", input.display()))?;
        results.push(stats);
    }

    if !args.common.json {
        emit_stdout(&results, args.tabular);
    }
    Ok(results)
}

fn emit_stdout(results: &[FastaStats], tabular: bool) {
    if tabular {
        for (i, s) in results.iter().enumerate() {
            let rendered = render_tabular(s);
            if i == 0 {
                print!("{rendered}");
            } else {
                for line in rendered.lines().skip(1) {
                    println!("{line}");
                }
            }
        }
    } else {
        for s in results {
            print!("{}", render_pretty(s));
        }
    }
}

fn main() -> ExitCode {
    let raw_args: Vec<String> = std::env::args().collect();
    if let Some(mode) = intercept_help(&raw_args) {
        render_help(&HELP, mode);
        return ExitCode::SUCCESS;
    }
    let args = Cli::parse();
    let common = args.common.clone();
    run(&common, META, || pipeline(&args))
}
