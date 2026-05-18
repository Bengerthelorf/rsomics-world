use rsomics_common::{Result, RsomicsError};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

pub fn sample(input: &Path, fraction: f64, seed: u64, output: &mut dyn Write) -> Result<u64> {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};
    let file = File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let reader = BufReader::new(file);
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    let mut count: u64 = 0;
    let mut idx: u64 = 0;
    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') {
            writeln!(out, "{line}").map_err(RsomicsError::Io)?;
            continue;
        }
        let mut h = DefaultHasher::new();
        (seed, idx).hash(&mut h);
        let hash_val = h.finish();
        let threshold = (fraction * u64::MAX as f64) as u64;
        if hash_val < threshold {
            writeln!(out, "{line}").map_err(RsomicsError::Io)?;
            count += 1;
        }
        idx += 1;
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
