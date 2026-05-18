use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub fn makewindows(
    genome_path: &Path,
    window_size: u64,
    step: u64,
    output: &mut dyn Write,
) -> Result<u64> {
    let file = File::open(genome_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", genome_path.display())))?;
    let reader = BufReader::new(file);
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    let mut count: u64 = 0;

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 2 {
            continue;
        }
        let chrom = f[0];
        let chrom_len: u64 = f[1]
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("bad length: {e}")))?;

        let mut start = 0u64;
        while start < chrom_len {
            let end = (start + window_size).min(chrom_len);
            writeln!(out, "{chrom}\t{start}\t{end}").map_err(RsomicsError::Io)?;
            count += 1;
            start += step;
        }
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
