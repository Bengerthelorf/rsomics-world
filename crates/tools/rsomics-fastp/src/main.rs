use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;

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
}

fn main() -> Result<()> {
    let args = Args::parse();
    rsomics_fastp::io::copy_se(&args.in1, &args.out1)
}
