use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub fn gff_to_bed(input: &Path, output: &mut dyn Write) -> Result<u64> {
    let file = File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let reader = BufReader::new(file);
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    let mut count: u64 = 0;

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        let chrom = fields[0];
        let start: u64 = fields[3]
            .parse::<u64>()
            .map_err(|e| RsomicsError::InvalidInput(format!("bad start: {e}")))?
            .saturating_sub(1); // GFF is 1-based, BED is 0-based
        let end: u64 = fields[4]
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("bad end: {e}")))?;
        let name = fields[2];
        let score = fields[5];
        let strand = fields[6];

        writeln!(out, "{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}")
            .map_err(RsomicsError::Io)?;
        count += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
