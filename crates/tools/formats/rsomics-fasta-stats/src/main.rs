mod cli;

use std::process::ExitCode;

use clap::Parser;
use rsomics_common::run;
use rsomics_help::{intercept_help, render as render_help};

use cli::{Cli, HELP, META};

fn main() -> ExitCode {
    let raw_args: Vec<String> = std::env::args().collect();
    if let Some(mode) = intercept_help(&raw_args) {
        render_help(&HELP, mode);
        return ExitCode::SUCCESS;
    }
    let args = Cli::parse();
    let common = args.common.clone();
    run(&common, META, || args.execute())
}
