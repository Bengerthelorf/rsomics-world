use std::process::ExitCode;

use clap::{Parser, Subcommand};
use rsomics_common::{CommonFlags, Result, StderrLog, run};

use rsomics_bam::cmd::view::{self, ViewArgs};

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-bam",
    version,
    about = "Pure-Rust CLI over rust-htslib for BAM operations. FFI-wrapper, output-compatible with samtools.",
    long_about = None,
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,

    #[command(flatten)]
    common: CommonFlags,
}

#[derive(Debug, Subcommand)]
enum Cmd {
    /// View / count BAM records (analogous to `samtools view`).
    View(ViewArgs),
}

fn dispatch(cli: &Cli) -> Result<()> {
    let log = StderrLog::from_flags(&cli.common);
    match &cli.cmd {
        Cmd::View(a) => view::run(a, &log),
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let common = cli.common.clone();
    run(&common, || dispatch(&cli))
}
