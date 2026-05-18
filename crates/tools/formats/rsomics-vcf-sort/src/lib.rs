use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub fn sort_vcf(input: &Path, output: &mut dyn Write) -> Result<u64> {
    let file = File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let reader = BufReader::new(file);
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    let mut header_lines = Vec::new();
    let mut records: Vec<(String, u64, String)> = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') {
            header_lines.push(line);
            continue;
        }
        let pos = line
            .split('\t')
            .nth(1)
            .and_then(|s| s.parse::<u64>().ok())
            .unwrap_or(0);
        let chrom = line.split('\t').next().unwrap_or("").to_string();
        records.push((chrom, pos, line));
    }

    records.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    for h in &header_lines {
        writeln!(out, "{h}").map_err(RsomicsError::Io)?;
    }
    for (_, _, line) in &records {
        writeln!(out, "{line}").map_err(RsomicsError::Io)?;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(records.len() as u64)
}
