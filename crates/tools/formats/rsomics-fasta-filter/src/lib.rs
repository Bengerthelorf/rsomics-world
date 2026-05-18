use std::io::{BufWriter, Write};
use std::path::Path;

use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};

pub fn filter(input: &Path, min_len: usize, max_len: usize, output: &mut dyn Write) -> Result<u64> {
    let mut reader = parse_fastx_file(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut count: u64 = 0;

    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| RsomicsError::InvalidInput(format!("parsing: {e}")))?;
        let len = rec.seq().len();
        if len >= min_len && (max_len == 0 || len <= max_len) {
            out.write_all(b">").map_err(RsomicsError::Io)?;
            out.write_all(rec.id()).map_err(RsomicsError::Io)?;
            out.write_all(b"\n").map_err(RsomicsError::Io)?;
            out.write_all(&rec.seq()).map_err(RsomicsError::Io)?;
            out.write_all(b"\n").map_err(RsomicsError::Io)?;
            count += 1;
        }
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
