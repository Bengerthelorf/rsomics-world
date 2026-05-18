use std::io::{BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::open_fastq;

fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'N' | b'n' => b'N',
        other => other,
    }
}

pub fn revcomp(input: &Path, output: &mut dyn Write) -> Result<u64> {
    let mut reader = open_fastq(input)?;
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut count: u64 = 0;

    for result in reader.by_ref() {
        let rec = result?;
        let rc_seq: Vec<u8> = rec.seq.iter().rev().map(|&b| complement(b)).collect();
        let rev_qual: Vec<u8> = rec.qual.iter().rev().copied().collect();

        out.write_all(b"@").map_err(RsomicsError::Io)?;
        out.write_all(&rec.id).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
        out.write_all(&rc_seq).map_err(RsomicsError::Io)?;
        out.write_all(b"\n+\n").map_err(RsomicsError::Io)?;
        out.write_all(&rev_qual).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
        count += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
