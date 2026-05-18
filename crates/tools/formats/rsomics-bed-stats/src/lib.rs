#![allow(clippy::cast_precision_loss)]

use rsomics_common::{Result, RsomicsError};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub struct BedStats {
    pub count: u64,
    pub total_bases: u64,
    pub min_len: u64,
    pub max_len: u64,
    pub mean_len: f64,
    pub median_len: f64,
}

pub fn stats(input: &Path) -> Result<BedStats> {
    let file = File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let reader = BufReader::new(file);
    let mut lengths: Vec<u64> = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        let s: u64 = f[1]
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("start: {e}")))?;
        let e: u64 = f[2]
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("end: {e}")))?;
        lengths.push(e.saturating_sub(s));
    }
    if lengths.is_empty() {
        return Ok(BedStats {
            count: 0,
            total_bases: 0,
            min_len: 0,
            max_len: 0,
            mean_len: 0.0,
            median_len: 0.0,
        });
    }
    lengths.sort_unstable();
    let total: u64 = lengths.iter().sum();
    let n = lengths.len();
    let median = if n % 2 == 0 {
        (lengths[n / 2 - 1] + lengths[n / 2]) as f64 / 2.0
    } else {
        lengths[n / 2] as f64
    };
    Ok(BedStats {
        count: n as u64,
        total_bases: total,
        min_len: lengths[0],
        max_len: lengths[n - 1],
        mean_len: total as f64 / n as f64,
        median_len: median,
    })
}
