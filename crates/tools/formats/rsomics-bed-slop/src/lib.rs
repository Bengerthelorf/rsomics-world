use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub fn slop(
    bed_path: &Path,
    genome_path: &Path,
    left: u64,
    right: u64,
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
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        let chrom = fields[0];
        let start: u64 = fields[1]
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("bad start: {e}")))?;
        let end: u64 = fields[2]
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("bad end: {e}")))?;

        let chrom_len = genome.get(chrom).copied().unwrap_or(u64::MAX);
        let new_start = start.saturating_sub(left);
        let new_end = (end + right).min(chrom_len);

        writeln!(out, "{chrom}\t{new_start}\t{new_end}").map_err(RsomicsError::Io)?;
        count += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}

fn load_genome(path: &Path) -> Result<HashMap<String, u64>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("genome file {}: {e}", path.display())))?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 2 {
            let len: u64 = fields[1]
                .parse()
                .map_err(|e| RsomicsError::InvalidInput(format!("bad genome length: {e}")))?;
            map.insert(fields[0].to_string(), len);
        }
    }

    Ok(map)
}
