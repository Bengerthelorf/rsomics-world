use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub fn concat(inputs: &[&Path], output: &mut dyn Write) -> Result<u64> {
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    let mut header_written = false;
    let mut total: u64 = 0;

    for input in inputs {
        let file = File::open(input)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line.map_err(RsomicsError::Io)?;
            if line.starts_with('#') {
                if !header_written {
                    writeln!(out, "{line}").map_err(RsomicsError::Io)?;
                }
                continue;
            }
            header_written = true;
            writeln!(out, "{line}").map_err(RsomicsError::Io)?;
            total += 1;
        }
        header_written = true;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(total)
}
