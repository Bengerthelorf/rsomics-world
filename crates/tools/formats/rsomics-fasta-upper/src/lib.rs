use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn uppercase(input: &Path, output: &mut dyn Write) -> Result<u64> {
    let mut reader = parse_fastx_file(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut count: u64 = 0;
    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| RsomicsError::InvalidInput(format!("parsing: {e}")))?;
        out.write_all(b">").map_err(RsomicsError::Io)?;
        out.write_all(rec.id()).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
        let upper: Vec<u8> = rec.seq().iter().map(u8::to_ascii_uppercase).collect();
        out.write_all(&upper).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
        count += 1;
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
