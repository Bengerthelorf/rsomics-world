use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use rsomics_common::{CommonFlags, Context, Result, RsomicsError, ToolMeta, run};
use rsomics_help::{HelpMode, intercept_help};

use rsomics_fasta_stats::{Config, FastaStats, compute_stats, render_pretty, render_tabular};

const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

const TAGLINE: &str =
    "Per-file statistics for FASTA inputs (Rust port of seqkit stats — FASTA subset).";

/// Compute per-file statistics for FASTA inputs.
///
/// Output is the FASTA subset of `seqkit stats`. For FASTQ-only quality
/// columns, use the sibling FASTQ-stats crate.
// `disable_help_flag` lets us intercept `--help` / `-h` from raw argv
// before clap parses, so `--help --plain` and `--help --json` can route
// to the renderer in `rsomics-help` instead of clap's classic surface.
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

fn print_rich_help() {
    use rsomics_help::{Banner, FlagRowSpec, example_line, flag_table, section_header, tagline};
    let color = !rsomics_help::no_color_env();
    println!();
    println!("{}", Banner::family(META.name).render(color));
    println!();
    println!("  {}", tagline(META.name, META.version, TAGLINE, color));
    println!();
    println!("{}", section_header("USAGE", color));
    println!("  rsomics-fasta-stats [OPTIONS] <INPUTS>...");
    println!();
    println!("{}", section_header("OPTIONS", color));
    println!(
        "{}",
        flag_table(
            &[
                FlagRowSpec {
                    short: Some('a'),
                    long: "all",
                    value: None,
                    desc: "Emit extended stats (Q1/Q2/Q3, N50, GC%, sum_gap, sum_n)",
                },
                FlagRowSpec {
                    short: Some('T'),
                    long: "tabular",
                    value: None,
                    desc: "Tab-separated machine-readable output",
                },
                FlagRowSpec {
                    short: Some('G'),
                    long: "gap-letters",
                    value: Some("<CHARS>"),
                    desc: "Gap chars for --all (default \"- .\", matches seqkit)",
                },
                FlagRowSpec {
                    short: None,
                    long: "json",
                    value: None,
                    desc: "Emit AI-friendly JSON envelope on stdout",
                },
                FlagRowSpec {
                    short: Some('t'),
                    long: "threads",
                    value: Some("<N>"),
                    desc: "Worker thread count (default: available cores)",
                },
                FlagRowSpec {
                    short: Some('h'),
                    long: "help",
                    value: None,
                    desc: "Show this help; add `--plain` or `--json` for alt modes",
                },
            ],
            color,
        )
    );
    println!();
    println!("{}", section_header("EXAMPLES", color));
    println!(
        "{}",
        example_line("Default stats", "rsomics-fasta-stats genome.fa", color)
    );
    println!(
        "{}",
        example_line(
            "Extended (--all), tabular, on a gzip input",
            "rsomics-fasta-stats --tabular --all genome.fa.gz",
            color
        )
    );
    println!(
        "{}",
        example_line(
            "JSON envelope piped through jq",
            "rsomics-fasta-stats --json scaffolds.fasta | jq .result",
            color
        )
    );
    println!();
}

fn print_json_help() {
    use rsomics_help::{Example, FlagGroup, FlagSpec, HelpJson, Origin};
    let help = HelpJson {
        origin: Some(Origin {
            upstream: "seqkit",
            upstream_license: "MIT",
            our_license: "MIT OR Apache-2.0",
            paper_doi: Some("10.1371/journal.pone.0163962"),
        }),
        flag_groups: vec![FlagGroup {
            title: "Options",
            flags: vec![
                FlagSpec {
                    short: Some('a'),
                    long: "all",
                    aliases: vec![],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "Emit extended statistics block",
                    why_default: Some(
                        "seqkit's `--all` opts in to the extended columns; off by default",
                    ),
                },
                FlagSpec {
                    short: Some('T'),
                    long: "tabular",
                    aliases: vec![],
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
                    aliases: vec![],
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
                    aliases: vec![],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "Emit AI-friendly JSON envelope on stdout",
                    why_default: None,
                },
            ],
        }],
        examples: vec![
            Example {
                description: "Default stats",
                command: "rsomics-fasta-stats genome.fa",
            },
            Example {
                description: "Extended tabular on gzip input",
                command: "rsomics-fasta-stats --tabular --all genome.fa.gz",
            },
        ],
        json_result_schema_doc: Some(
            "https://docs.rs/rsomics-fasta-stats/0.2.0/#json-output-schema",
        ),
        ..HelpJson::new(META.name, META.version, TAGLINE)
    };
    let _ = serde_json::to_writer_pretty(std::io::stdout().lock(), &help);
    println!();
}

fn print_plain_help() {
    use rsomics_help::{FlagRowSpec, flag_table};
    println!("{} {} — {}", META.name, META.version, TAGLINE);
    println!();
    println!("USAGE");
    println!("  rsomics-fasta-stats [OPTIONS] <INPUTS>...");
    println!();
    println!("OPTIONS");
    println!(
        "{}",
        flag_table(
            &[
                FlagRowSpec {
                    short: Some('a'),
                    long: "all",
                    value: None,
                    desc: "Emit extended stats (Q1/Q2/Q3, N50, GC%, sum_gap, sum_n)",
                },
                FlagRowSpec {
                    short: Some('T'),
                    long: "tabular",
                    value: None,
                    desc: "Tab-separated machine-readable output",
                },
                FlagRowSpec {
                    short: Some('G'),
                    long: "gap-letters",
                    value: Some("<CHARS>"),
                    desc: "Gap chars for --all (default \"- .\")",
                },
                FlagRowSpec {
                    short: None,
                    long: "json",
                    value: None,
                    desc: "Emit AI-friendly JSON envelope on stdout",
                },
                FlagRowSpec {
                    short: Some('t'),
                    long: "threads",
                    value: Some("<N>"),
                    desc: "Worker thread count",
                },
                FlagRowSpec {
                    short: Some('h'),
                    long: "help",
                    value: None,
                    desc: "Show this help (--help for rich, --help --json for AI)",
                },
            ],
            false,
        )
    );
    println!();
    println!("EXAMPLES");
    println!("  rsomics-fasta-stats genome.fa                       # Default stats");
    println!("  rsomics-fasta-stats --tabular --all genome.fa.gz    # Extended, tabular");
    println!("  rsomics-fasta-stats --json scaffolds.fa | jq .result  # JSON pipe");
    println!();
}

fn main() -> ExitCode {
    let raw_args: Vec<String> = std::env::args().collect();
    if let Some(mode) = intercept_help(&raw_args) {
        match mode {
            HelpMode::Rich => print_rich_help(),
            HelpMode::Plain => print_plain_help(),
            HelpMode::Json => print_json_help(),
        }
        return ExitCode::SUCCESS;
    }
    let args = Cli::parse();
    let common = args.common.clone();
    run(&common, META, || pipeline(&args))
}
