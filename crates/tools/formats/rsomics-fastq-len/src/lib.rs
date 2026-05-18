use std::io::{BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::open_fastq;

pub fn lengths(input: &Path, tab: bool, output: &mut dyn Write) -> Result<u64> {
    let mut reader = open_fastq(input)?;
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut count: u64 = 0;

    for result in reader.by_ref() {
        let rec = result?;
        let len = rec.seq.len();
        if tab {
            out.write_all(&rec.id).map_err(RsomicsError::Io)?;
            writeln!(out, "\t{len}").map_err(RsomicsError::Io)?;
        } else {
            writeln!(out, "{len}").map_err(RsomicsError::Io)?;
        }
        count += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
