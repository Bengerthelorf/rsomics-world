#![allow(clippy::cast_precision_loss)]
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

struct Interval {
    chrom: String,
    start: u64,
    end: u64,
}

fn load_bed(path: &Path) -> Result<Vec<Interval>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let reader = BufReader::new(file);
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        out.push(Interval {
            chrom: f[0].to_string(),
            start: f[1]
                .parse()
                .map_err(|e| RsomicsError::InvalidInput(format!("start: {e}")))?,
            end: f[2]
                .parse()
                .map_err(|e| RsomicsError::InvalidInput(format!("end: {e}")))?,
        });
    }
    out.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));
    Ok(out)
}

#[derive(Debug)]
pub struct JaccardResult {
    pub intersection: u64,
    pub union: u64,
    pub jaccard: f64,
    pub n_intersections: u64,
}

pub fn jaccard(a_path: &Path, b_path: &Path) -> Result<JaccardResult> {
    let a = load_bed(a_path)?;
    let b = load_bed(b_path)?;

    let mut intersection: u64 = 0;
    let mut n_intersections: u64 = 0;
    let mut j = 0usize;

    for ai in &a {
        while j < b.len()
            && (b[j].chrom < ai.chrom || (b[j].chrom == ai.chrom && b[j].end <= ai.start))
        {
            j += 1;
        }
        let mut k = j;
        while k < b.len() && b[k].chrom == ai.chrom && b[k].start < ai.end {
            let ov_start = ai.start.max(b[k].start);
            let ov_end = ai.end.min(b[k].end);
            if ov_start < ov_end {
                intersection += ov_end - ov_start;
                n_intersections += 1;
            }
            k += 1;
        }
    }

    let sum_a: u64 = a.iter().map(|i| i.end - i.start).sum();
    let sum_b: u64 = b.iter().map(|i| i.end - i.start).sum();
    let union = sum_a + sum_b - intersection;

    let jaccard_val = if union == 0 {
        0.0
    } else {
        intersection as f64 / union as f64
    };

    Ok(JaccardResult {
        intersection,
        union,
        jaccard: jaccard_val,
        n_intersections,
    })
}
