use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use rsomics_common::{CommonFlags, Context, Result, RsomicsError, ToolMeta, run};

use rsomics_fasta_stats::{Config, FastaStats, compute_stats, render_pretty, render_tabular};

const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

/// Compute per-file statistics for FASTA inputs.
///
/// Output is the FASTA subset of `seqkit stats`. For FASTQ-only quality
/// columns, use the sibling FASTQ-stats crate.
#[derive(Parser, Debug)]
#[command(name = "rsomics-fasta-stats", version, about, long_about = None)]
struct Cli {
    /// FASTA file(s). Use `-` for stdin. Gzip, bzip2, xz, zstd inputs are
    /// auto-detected by extension or magic bytes.
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
    let args = Cli::parse();
    let common = args.common.clone();
    run(&common, META, || pipeline(&args))
}
