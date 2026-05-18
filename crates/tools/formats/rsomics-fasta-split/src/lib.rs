use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};

pub fn split_by_count(input: &Path, seqs_per_file: u64, out_prefix: &str) -> Result<u64> {
    let mut reader = parse_fastx_file(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;

    let mut file_idx: u64 = 0;
    let mut count_in_file: u64 = 0;
    let mut writer: Option<BufWriter<File>> = None;

    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| RsomicsError::InvalidInput(format!("parsing: {e}")))?;

        if count_in_file >= seqs_per_file || writer.is_none() {
            let path = format!("{out_prefix}{file_idx:04}.fa");
            let f = File::create(&path)
                .map_err(|e| RsomicsError::InvalidInput(format!("creating {path}: {e}")))?;
            writer = Some(BufWriter::new(f));
            file_idx += 1;
            count_in_file = 0;
        }

        let w = writer.as_mut().unwrap();
        w.write_all(b">").map_err(RsomicsError::Io)?;
        w.write_all(rec.id()).map_err(RsomicsError::Io)?;
        w.write_all(b"\n").map_err(RsomicsError::Io)?;
        w.write_all(&rec.seq()).map_err(RsomicsError::Io)?;
        w.write_all(b"\n").map_err(RsomicsError::Io)?;
        count_in_file += 1;
    }

    if let Some(mut w) = writer {
        w.flush().map_err(RsomicsError::Io)?;
    }
    Ok(file_idx)
}
