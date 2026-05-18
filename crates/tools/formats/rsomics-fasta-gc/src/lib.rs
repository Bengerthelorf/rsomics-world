#![allow(clippy::cast_precision_loss)]
use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn gc_content(input: &Path, output: &mut dyn Write) -> Result<u64> {
    let mut reader = parse_fastx_file(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    let mut count: u64 = 0;
    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| RsomicsError::InvalidInput(format!("parsing: {e}")))?;
        let seq = rec.seq();
        let gc = seq
            .iter()
            .filter(|b| matches!(b, b'G' | b'C' | b'g' | b'c'))
            .count();
        let pct = if seq.is_empty() {
            0.0
        } else {
            gc as f64 / seq.len() as f64 * 100.0
        };
        out.write_all(rec.id()).map_err(RsomicsError::Io)?;
        writeln!(out, "\t{:.2}", pct).map_err(RsomicsError::Io)?;
        count += 1;
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
