use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use needletail::parse_fastx_file;

/// Identity-copy a single-end FASTQ file. No transformation; validates the
/// reader / writer plumbing in isolation before filtering layers ride on top.
///
/// # Errors
///
/// Returns `Err` if the input cannot be opened, a record fails to parse,
/// the output cannot be created, or a write to the output fails.
pub fn copy_se(input: &Path, output: &Path) -> Result<()> {
    let mut reader = parse_fastx_file(input)
        .with_context(|| format!("opening input FASTQ {}", input.display()))?;
    let mut writer = BufWriter::new(
        File::create(output)
            .with_context(|| format!("creating output FASTQ {}", output.display()))?,
    );
    while let Some(record) = reader.next() {
        let rec = record.context("malformed FASTQ record")?;
        rec.write(&mut writer, None)
            .context("writing record to output")?;
    }
    writer.flush().context("flushing output writer")?;
    Ok(())
}
