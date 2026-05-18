use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use regex::Regex;
use rsomics_common::{Result, RsomicsError};

pub fn grep(input: &Path, pattern: &str, invert: bool, output: &mut dyn Write) -> Result<u64> {
    let re =
        Regex::new(pattern).map_err(|e| RsomicsError::InvalidInput(format!("bad regex: {e}")))?;
    let file = File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let reader = BufReader::new(file);
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    let mut count: u64 = 0;

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') {
            writeln!(out, "{line}").map_err(RsomicsError::Io)?;
            continue;
        }
        let matches = re.is_match(&line);
        if matches != invert {
            writeln!(out, "{line}").map_err(RsomicsError::Io)?;
            count += 1;
        }
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
