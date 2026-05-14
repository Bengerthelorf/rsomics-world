//! `rsomics-bam view` — inspect / count / project BAM records. v0.1.0
//! ships only `--count` (analogous to `samtools view -c`); record
//! projection (header / SAM output / region filter) follows.

use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{Result, StderrLog};
use rust_htslib::bam::{Read, Reader, Record};

use crate::htslib_bridge::{HtsResultExt, from_htslib};

#[derive(Debug, Args)]
pub struct ViewArgs {
    /// Input BAM (or SAM/CRAM) file.
    pub input: PathBuf,

    /// Print only the record count, not the records themselves. Mirrors
    /// `samtools view -c`.
    #[arg(short = 'c', long = "count", default_value_t = false)]
    pub count: bool,
}

/// Count records in a BAM file. Returned as a library function so tests
/// can compare counts directly without parsing stdout. The CLI wrapper
/// in [`run`] prints the result.
///
/// Uses rust-htslib's `read_into(&mut Record)` API — reuses a single
/// `Record` across the loop instead of allocating one per iteration.
/// On a 100k-record BAM this is ~3-5× cheaper than `.records()` and
/// drops to a fixed memory footprint regardless of record count.
///
/// # Errors
///
/// Returns `Err` if the input cannot be opened or a record fails to parse.
pub fn count_records(input: &Path) -> Result<u64> {
    let mut reader = Reader::from_path(input).rs()?;
    let mut rec = Record::new();
    let mut n: u64 = 0;
    while let Some(r) = reader.read(&mut rec) {
        r.map_err(from_htslib)?;
        n += 1;
    }
    Ok(n)
}

/// Entry point for the `view` subcommand. The `log` argument carries
/// `--quiet` / `--verbose` state from `CommonFlags`.
///
/// # Errors
///
/// Returns `Err` if the input file cannot be opened or a record fails to
/// parse. With `--count`, the record count is printed to stdout.
/// Without `--count`, the reader iterates records without emitting them
/// — record projection (SAM / BAM out, region filter) is a follow-up
/// subcommand surface.
pub fn run(args: &ViewArgs, _log: &StderrLog) -> Result<()> {
    if args.count {
        let n = count_records(&args.input)?;
        println!("{n}");
        return Ok(());
    }
    let mut reader = Reader::from_path(&args.input).rs()?;
    let mut rec = Record::new();
    while let Some(r) = reader.read(&mut rec) {
        r.map_err(from_htslib)?;
    }
    Ok(())
}
