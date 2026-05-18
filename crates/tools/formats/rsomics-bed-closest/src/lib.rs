use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

#[derive(Debug, Clone)]
struct Interval {
    chrom: String,
    start: u64,
    end: u64,
}

fn load_sorted_bed(path: &Path) -> Result<Vec<Interval>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let reader = BufReader::new(file);
    let mut intervals = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        intervals.push(Interval {
            chrom: fields[0].to_string(),
            start: fields[1]
                .parse()
                .map_err(|e| RsomicsError::InvalidInput(format!("bad start: {e}")))?,
            end: fields[2]
                .parse()
                .map_err(|e| RsomicsError::InvalidInput(format!("bad end: {e}")))?,
        });
    }
    Ok(intervals)
}

pub fn closest(a_path: &Path, b_path: &Path, output: &mut dyn Write) -> Result<u64> {
    let a_intervals = load_sorted_bed(a_path)?;
    let b_intervals = load_sorted_bed(b_path)?;
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    let mut count: u64 = 0;

    for a in &a_intervals {
        let mut best_dist: i64 = i64::MAX;
        let mut best_b: Option<&Interval> = None;

        for b in &b_intervals {
            if a.chrom != b.chrom {
                continue;
            }
            let dist = if a.end <= b.start {
                (b.start - a.end) as i64
            } else if b.end <= a.start {
                (a.start - b.end) as i64
            } else {
                0
            };
            if dist < best_dist {
                best_dist = dist;
                best_b = Some(b);
            }
        }

        if let Some(b) = best_b {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{best_dist}",
                a.chrom, a.start, a.end, b.chrom, b.start, b.end
            )
            .map_err(RsomicsError::Io)?;
        } else {
            writeln!(out, "{}\t{}\t{}\t.\t-1\t-1\t-1", a.chrom, a.start, a.end)
                .map_err(RsomicsError::Io)?;
        }
        count += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
