use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use rsomics_fasta_index::fetch_region;

pub fn getfasta(bed_path: &Path, fasta_path: &Path, output: &mut dyn Write) -> Result<u64> {
    let fai_path = fasta_path.with_extension(format!(
        "{}.fai",
        fasta_path.extension().unwrap_or_default().to_string_lossy()
    ));

    let file = File::open(bed_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bed_path.display())))?;
    let reader = BufReader::new(file);
    let mut out = BufWriter::with_capacity(256 * 1024, output);
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
        let start = fields[1];
        let end = fields[2];
        let region = format!("{chrom}:{start}-{end}");

        let seq = fetch_region(fasta_path, &fai_path, &region)?;
        writeln!(out, ">{region}").map_err(RsomicsError::Io)?;
        out.write_all(&seq).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
        count += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
