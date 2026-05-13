use anyhow::Result;
use clap::Parser;

#[derive(Debug, Parser)]
#[command(name = "rsomics-fastp", version, about, long_about = None)]
struct Args {}

fn main() -> Result<()> {
    let _args = Args::parse();
    anyhow::bail!("rsomics-fastp is in early scaffold; no preprocessing implemented yet")
}
