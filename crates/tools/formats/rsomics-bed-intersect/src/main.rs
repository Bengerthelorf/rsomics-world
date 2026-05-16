mod cli;

use std::process::{self, ExitCode};

use clap::Parser;
use rsomics_common::{ExitCode as RsomicsExit, run};
use rsomics_help::{intercept_help, render as render_help};

use cli::{Cli, HELP, META};

fn main() -> ExitCode {
    let raw_args: Vec<String> = std::env::args().collect();
    if let Some(mode) = intercept_help(&raw_args) {
        render_help(&HELP, mode);
        return process::ExitCode::from(RsomicsExit::Ok);
    }
    let cli = Cli::parse();
    let common = cli.common.clone();
    run(&common, META, || cli.execute())
}
