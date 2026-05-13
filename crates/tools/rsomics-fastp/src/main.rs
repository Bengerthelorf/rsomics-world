use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;

use rsomics_fastp::filter::FilterConfig;

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-fastp",
    version,
    about = "Pure-Rust FASTQ preprocessor — output-compatible with OpenGene/fastp.",
    long_about = None,
)]
struct Args {
    /// Input FASTQ (single-end). Compression auto-detected by file content.
    #[arg(short = 'i', long = "in1")]
    in1: PathBuf,

    /// Output FASTQ.
    #[arg(short = 'o', long = "out1")]
    out1: PathBuf,

    /// JSON report path (fastp-compatible schema subset).
    #[arg(short = 'j', long = "json")]
    json: Option<PathBuf>,

    /// Per-base quality threshold (Phred); bases below this count as unqualified.
    #[arg(long = "qualified_quality_phred", default_value_t = 15)]
    qualified_quality_phred: u8,

    /// Reject a read if the fraction of unqualified bases (percent) exceeds this.
    #[arg(long = "unqualified_percent_limit", default_value_t = 40)]
    unqualified_percent_limit: u32,

    /// Reject reads shorter than this many bases.
    #[arg(long = "length_required", default_value_t = 15)]
    length_required: usize,

    /// Reject reads with more than this many N bases.
    #[arg(long = "n_base_limit", default_value_t = 5)]
    n_base_limit: usize,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let cfg = FilterConfig {
        qualified_quality_phred: args.qualified_quality_phred,
        unqualified_percent_limit: args.unqualified_percent_limit,
        length_required: args.length_required,
        n_base_limit: args.n_base_limit,
    };
    rsomics_fastp::io::process_se(&args.in1, &args.out1, args.json.as_deref(), cfg)?;
    Ok(())
}
