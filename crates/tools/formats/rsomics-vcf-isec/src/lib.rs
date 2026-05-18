use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

fn variant_key(line: &str) -> Option<String> {
    let fields: Vec<&str> = line.splitn(6, '\t').collect();
    if fields.len() < 5 {
        return None;
    }
    Some(format!(
        "{}\t{}\t{}\t{}",
        fields[0], fields[1], fields[3], fields[4]
    ))
}

fn load_keys(path: &Path) -> Result<HashSet<String>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let reader = BufReader::new(file);
    let mut keys = HashSet::new();
    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') {
            continue;
        }
        if let Some(k) = variant_key(&line) {
            keys.insert(k);
        }
    }
    Ok(keys)
}

pub fn isec(a_path: &Path, b_path: &Path, output: &mut dyn Write) -> Result<u64> {
    let b_keys = load_keys(b_path)?;
    let file = File::open(a_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", a_path.display())))?;
    let reader = BufReader::new(file);
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    let mut count: u64 = 0;

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') {
            writeln!(out, "{line}").map_err(RsomicsError::Io)?;
            continue;
        }
        if let Some(k) = variant_key(&line) {
            if b_keys.contains(&k) {
                writeln!(out, "{line}").map_err(RsomicsError::Io)?;
                count += 1;
            }
        }
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
