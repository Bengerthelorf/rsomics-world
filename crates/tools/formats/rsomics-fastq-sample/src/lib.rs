use std::io::{BufWriter, Write};
use std::path::Path;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::open_fastq;

pub fn sample(input: &Path, fraction: f64, seed: u64, output: &mut dyn Write) -> Result<u64> {
    let mut reader = open_fastq(input)?;
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut rng = StdRng::seed_from_u64(seed);
    let mut count: u64 = 0;

    for result in reader.by_ref() {
        let rec = result?;
        if rng.r#gen::<f64>() < fraction {
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
