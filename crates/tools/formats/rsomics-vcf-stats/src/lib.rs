#![allow(clippy::cast_precision_loss)]

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

#[derive(Debug, Default, Serialize)]
pub struct VcfStats {
    pub total: u64,
    pub snps: u64,
    pub insertions: u64,
    pub deletions: u64,
    pub mnps: u64,
    pub transitions: u64,
    pub transversions: u64,
    pub ti_tv: f64,
}

fn is_transition(r: u8, a: u8) -> bool {
    matches!(
        (r, a),
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')
    )
}

pub fn stats(input: &Path) -> Result<VcfStats> {
    let file = File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let reader = BufReader::new(file);
    let mut s = VcfStats::default();

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.splitn(6, '\t').collect();
        if fields.len() < 5 {
            continue;
        }
        let ref_allele = fields[3];
        let alt_field = fields[4];

        for alt in alt_field.split(',') {
            if alt == "." || alt == "*" {
                continue;
            }
            s.total += 1;
            let rlen = ref_allele.len();
            let alen = alt.len();

            if rlen == 1 && alen == 1 {
                s.snps += 1;
                let rb = ref_allele.as_bytes()[0].to_ascii_uppercase();
                let ab = alt.as_bytes()[0].to_ascii_uppercase();
                if is_transition(rb, ab) {
                    s.transitions += 1;
                } else {
                    s.transversions += 1;
                }
            } else if rlen > alen {
                s.deletions += 1;
            } else if alen > rlen {
                s.insertions += 1;
            } else {
                s.mnps += 1;
            }
        }
    }

    s.ti_tv = if s.transversions > 0 {
        s.transitions as f64 / s.transversions as f64
    } else {
        0.0
    };

    Ok(s)
}
