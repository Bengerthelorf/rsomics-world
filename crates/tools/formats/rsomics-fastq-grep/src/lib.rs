use std::io::{BufWriter, Write};
use std::path::Path;

use regex::bytes::Regex;
use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::open_fastq;

pub fn grep(input: &Path, pattern: &str, invert: bool, output: &mut dyn Write) -> Result<u64> {
    let re =
        Regex::new(pattern).map_err(|e| RsomicsError::InvalidInput(format!("bad regex: {e}")))?;
    let mut reader = open_fastq(input)?;
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut count: u64 = 0;

    for result in reader.by_ref() {
        let rec = result?;
        let matches = re.is_match(&rec.id);
        if matches != invert {
            out.write_all(b"@").map_err(RsomicsError::Io)?;
            out.write_all(&rec.id).map_err(RsomicsError::Io)?;
            out.write_all(b"\n").map_err(RsomicsError::Io)?;
            out.write_all(&rec.seq).map_err(RsomicsError::Io)?;
            out.write_all(b"\n+\n").map_err(RsomicsError::Io)?;
            out.write_all(&rec.qual).map_err(RsomicsError::Io)?;
            out.write_all(b"\n").map_err(RsomicsError::Io)?;
            count += 1;
        }
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
