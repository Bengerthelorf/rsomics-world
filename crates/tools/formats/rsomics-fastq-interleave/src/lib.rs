use std::io::{BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::open_fastq;

fn write_rec(out: &mut impl Write, id: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
    out.write_all(b"@").map_err(RsomicsError::Io)?;
    out.write_all(id).map_err(RsomicsError::Io)?;
    out.write_all(b"\n").map_err(RsomicsError::Io)?;
    out.write_all(seq).map_err(RsomicsError::Io)?;
    out.write_all(b"\n+\n").map_err(RsomicsError::Io)?;
    out.write_all(qual).map_err(RsomicsError::Io)?;
    out.write_all(b"\n").map_err(RsomicsError::Io)
}

pub fn interleave(r1_path: &Path, r2_path: &Path, output: &mut dyn Write) -> Result<u64> {
    let mut r1 = open_fastq(r1_path)?;
    let mut r2 = open_fastq(r2_path)?;
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut pairs: u64 = 0;

    loop {
        match (r1.next(), r2.next()) {
            (Some(a), Some(b)) => {
                let a = a?;
                let b = b?;
                write_rec(&mut out, &a.id, &a.seq, &a.qual)?;
                write_rec(&mut out, &b.id, &b.seq, &b.qual)?;
                pairs += 1;
            }
            (None, None) => break,
            _ => {
                return Err(RsomicsError::InvalidInput(
                    "R1 and R2 have different record counts".into(),
                ));
            }
        }
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(pairs)
}
