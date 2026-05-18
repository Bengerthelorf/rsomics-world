use std::io::{BufWriter, Write};
use std::path::Path;

use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};

#[derive(Debug, Clone, Copy)]
pub enum SortKey {
    Name,
    Length,
    LengthDesc,
}

pub fn sort(input: &Path, key: SortKey, output: &mut dyn Write) -> Result<u64> {
    let mut reader = parse_fastx_file(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;

    let mut records: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| RsomicsError::InvalidInput(format!("parsing: {e}")))?;
        records.push((rec.id().to_vec(), rec.seq().to_vec()));
    }

    match key {
        SortKey::Name => records.sort_by(|a, b| a.0.cmp(&b.0)),
        SortKey::Length => records.sort_by_key(|r| r.1.len()),
        SortKey::LengthDesc => records.sort_by_key(|r| std::cmp::Reverse(r.1.len())),
    }

    let mut out = BufWriter::with_capacity(256 * 1024, output);
    for (id, seq) in &records {
        out.write_all(b">").map_err(RsomicsError::Io)?;
        out.write_all(id).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
        out.write_all(seq).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(records.len() as u64)
}
