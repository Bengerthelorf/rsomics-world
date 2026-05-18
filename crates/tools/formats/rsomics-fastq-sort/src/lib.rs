use std::io::{BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::{OwnedRecord, open_fastq};

#[derive(Debug, Clone, Copy)]
pub enum SortKey {
    Name,
    Length,
    LengthDesc,
}

pub fn sort(input: &Path, key: SortKey, output: &mut dyn Write) -> Result<u64> {
    let mut reader = open_fastq(input)?;
    let mut records: Vec<OwnedRecord> = Vec::new();

    for result in reader.by_ref() {
        records.push(result?);
    }

    match key {
        SortKey::Name => records.sort_by(|a, b| a.id.cmp(&b.id)),
        SortKey::Length => records.sort_by_key(|r| r.seq.len()),
        SortKey::LengthDesc => records.sort_by_key(|r| std::cmp::Reverse(r.seq.len())),
    }

    let mut out = BufWriter::with_capacity(256 * 1024, output);
    for rec in &records {
        out.write_all(b"@").map_err(RsomicsError::Io)?;
        out.write_all(&rec.id).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
        out.write_all(&rec.seq).map_err(RsomicsError::Io)?;
        out.write_all(b"\n+\n").map_err(RsomicsError::Io)?;
        out.write_all(&rec.qual).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(records.len() as u64)
}
