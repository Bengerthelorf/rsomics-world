use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};
use std::path::Path;

pub fn count(input: &Path) -> Result<u64> {
    let mut reader = parse_fastx_file(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let mut n: u64 = 0;
    while let Some(record) = reader.next() {
        record.map_err(|e| RsomicsError::InvalidInput(format!("parsing: {e}")))?;
        n += 1;
    }
    Ok(n)
}
