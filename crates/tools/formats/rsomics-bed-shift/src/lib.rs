use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub fn shift(
    bed_path: &Path,
    genome_path: &Path,
    offset: i64,
    output: &mut dyn Write,
) -> Result<u64> {
    let genome = load_genome(genome_path)?;
    let file = File::open(bed_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bed_path.display())))?;
    let reader = BufReader::new(file);
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    let mut count: u64 = 0;

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        let chrom = f[0];
        let start: i64 = f[1]
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("start: {e}")))?;
        let end: i64 = f[2]
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("end: {e}")))?;
        let chrom_len = genome.get(chrom).copied().unwrap_or(i64::MAX as u64) as i64;

        let new_start = (start + offset).max(0).min(chrom_len);
        let new_end = (end + offset).max(0).min(chrom_len);

        if new_start < new_end {
            writeln!(out, "{chrom}\t{new_start}\t{new_end}").map_err(RsomicsError::Io)?;
            count += 1;
        }
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}

fn load_genome(path: &Path) -> Result<HashMap<String, u64>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();
    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() >= 2 {
            let len: u64 = f[1]
                .parse()
                .map_err(|e| RsomicsError::InvalidInput(format!("len: {e}")))?;
            map.insert(f[0].to_string(), len);
        }
    }
    Ok(map)
}
