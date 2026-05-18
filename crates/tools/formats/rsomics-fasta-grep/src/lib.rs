use std::io::{BufWriter, Write};
use std::path::Path;

use needletail::parse_fastx_file;
use regex::bytes::Regex;
use rsomics_common::{Result, RsomicsError};

pub fn grep(input: &Path, pattern: &str, invert: bool, output: &mut dyn Write) -> Result<u64> {
    let re =
        Regex::new(pattern).map_err(|e| RsomicsError::InvalidInput(format!("bad regex: {e}")))?;
    let mut reader = parse_fastx_file(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut count: u64 = 0;

    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| RsomicsError::InvalidInput(format!("parsing: {e}")))?;
        let matches = re.is_match(rec.id());
        if matches != invert {
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
