use std::io::{BufWriter, Write};
use std::path::Path;

use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};

pub fn lengths(input: &Path, tab: bool, output: &mut dyn Write) -> Result<u64> {
    let mut reader = parse_fastx_file(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut count: u64 = 0;

    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| RsomicsError::InvalidInput(format!("parsing: {e}")))?;
        let len = rec.seq().len();
        if tab {
            out.write_all(rec.id()).map_err(RsomicsError::Io)?;
            writeln!(out, "\t{len}").map_err(RsomicsError::Io)?;
        } else {
            writeln!(out, "{len}").map_err(RsomicsError::Io)?;
        }
        count += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
