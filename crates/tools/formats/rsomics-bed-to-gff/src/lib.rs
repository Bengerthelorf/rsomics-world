use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub fn bed_to_gff(
    input: &Path,
    source: &str,
    feature_type: &str,
    output: &mut dyn Write,
) -> Result<u64> {
    let file = File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let reader = BufReader::new(file);
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    writeln!(out, "##gff-version 3").map_err(RsomicsError::Io)?;
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
        let start: u64 = f[1]
            .parse::<u64>()
            .map_err(|e| RsomicsError::InvalidInput(format!("start: {e}")))?
            + 1; // BED 0-based → GFF 1-based
        let end = f[2];
        let name = if f.len() > 3 { f[3] } else { "." };
        let score = if f.len() > 4 { f[4] } else { "." };
        let strand = if f.len() > 5 { f[5] } else { "." };

        writeln!(
            out,
            "{chrom}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t.\tName={name}"
        )
        .map_err(RsomicsError::Io)?;
        count += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
