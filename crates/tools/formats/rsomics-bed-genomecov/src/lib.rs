#![allow(clippy::cast_precision_loss)]
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub fn genomecov(bed_path: &Path, genome_path: &Path, output: &mut dyn Write) -> Result<()> {
    let genome = load_genome(genome_path)?;
    let file = File::open(bed_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bed_path.display())))?;
    let reader = BufReader::new(file);
    let mut out = BufWriter::with_capacity(64 * 1024, output);

    let mut depth_hist: HashMap<String, HashMap<u32, u64>> = HashMap::new();

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
        let start: u64 = f[1]
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("start: {e}")))?;
        let end: u64 = f[2]
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("end: {e}")))?;
        let len = end.saturating_sub(start);
        let hist = depth_hist.entry(chrom.to_string()).or_default();
        *hist.entry(1).or_insert(0) += len;
    }

    let mut chroms: Vec<&String> = genome.keys().collect();
    chroms.sort();

    for chrom in &chroms {
        let chrom_len = genome[*chrom];
        let covered = depth_hist
            .get(*chrom)
            .and_then(|h| h.get(&1))
            .copied()
            .unwrap_or(0);
        let uncovered = chrom_len.saturating_sub(covered);

        if uncovered > 0 {
            writeln!(
                out,
                "{chrom}\t0\t{uncovered}\t{chrom_len}\t{:.6}",
                uncovered as f64 / chrom_len as f64
            )
            .map_err(RsomicsError::Io)?;
        }
        if covered > 0 {
            writeln!(
                out,
                "{chrom}\t1\t{covered}\t{chrom_len}\t{:.6}",
                covered as f64 / chrom_len as f64
            )
            .map_err(RsomicsError::Io)?;
        }
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(())
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
