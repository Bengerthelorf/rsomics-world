use std::io::{BufWriter, Write};
use std::path::Path;

use needletail::parse_fastx_file;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rsomics_common::{Result, RsomicsError};

pub fn sample(input: &Path, fraction: f64, seed: u64, output: &mut dyn Write) -> Result<u64> {
    let mut reader = parse_fastx_file(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut rng = StdRng::seed_from_u64(seed);
    let mut count: u64 = 0;

    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| RsomicsError::InvalidInput(format!("parsing: {e}")))?;
        if rng.r#gen::<f64>() < fraction {
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
