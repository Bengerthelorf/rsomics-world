use std::process::ExitCode;

use clap::{Parser, Subcommand};
use rsomics_common::{CommonFlags, Result, StderrLog, ToolMeta, run};

use rsomics_bam::cmd::view::{self, ViewArgs, ViewSummary};

const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

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

/// `--json` envelope `result` body. The enum carries the subcommand
/// tag in its `subcommand` field so callers branch cleanly on it.
#[derive(Debug, serde::Serialize)]
#[serde(untagged)]
enum DispatchSummary {
    View(ViewSummary),
}

fn dispatch(cli: &Cli) -> Result<DispatchSummary> {
    let log = StderrLog::from_flags(&cli.common);
    match &cli.cmd {
        Cmd::View(a) => view::run(a, &log).map(DispatchSummary::View),
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let common = cli.common.clone();
    run(&common, META, || dispatch(&cli))
}
