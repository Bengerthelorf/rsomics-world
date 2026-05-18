use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub struct GffStats {
    pub total: u64,
    pub by_type: BTreeMap<String, u64>,
    pub by_source: BTreeMap<String, u64>,
    pub by_chrom: BTreeMap<String, u64>,
}

pub fn stats(input: &Path) -> Result<GffStats> {
    let file = File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let reader = BufReader::new(file);
    let mut s = GffStats {
        total: 0,
        by_type: BTreeMap::new(),
        by_source: BTreeMap::new(),
        by_chrom: BTreeMap::new(),
    };

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        s.total += 1;
        *s.by_chrom.entry(fields[0].to_string()).or_insert(0) += 1;
        *s.by_source.entry(fields[1].to_string()).or_insert(0) += 1;
        *s.by_type.entry(fields[2].to_string()).or_insert(0) += 1;
    }

    Ok(s)
}
