mod cli;
use clap::Parser;
use cli::{Cli, HELP, META};
use rsomics_common::run as runner;
use rsomics_help::{intercept_help, render as render_help};
use std::process::ExitCode;
fn main() -> ExitCode {
    let raw_args: Vec<String> = std::env::args().collect();
    if let Some(mode) = intercept_help(&raw_args) {
        render_help(&HELP, mode);
        return ExitCode::SUCCESS;
    }
    let args = Cli::parse();
    let common = args.common.clone();
    runner(&common, META, || args.execute())
}
