use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use regex::Regex;
use rsomics_common::{Result, RsomicsError};

pub fn filter_gff(
    input: &Path,
    feature_type: Option<&str>,
    pattern: Option<&str>,
    output: &mut dyn Write,
) -> Result<u64> {
    let re = pattern
        .map(|p| Regex::new(p).map_err(|e| RsomicsError::InvalidInput(format!("bad regex: {e}"))))
        .transpose()?;

    let file = File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let reader = BufReader::new(file);
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    let mut count: u64 = 0;

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') || line.is_empty() {
            writeln!(out, "{line}").map_err(RsomicsError::Io)?;
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        if feature_type.is_some_and(|ft| fields[2] != ft) {
            continue;
        }

        if re.as_ref().is_some_and(|r| !r.is_match(&line)) {
            continue;
        }

        writeln!(out, "{line}").map_err(RsomicsError::Io)?;
        count += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
