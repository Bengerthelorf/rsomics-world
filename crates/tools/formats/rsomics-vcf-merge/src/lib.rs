use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub fn merge(inputs: &[&Path], output: &mut dyn Write) -> Result<u64> {
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    let mut all_records: BTreeMap<(String, u64, String, String), Vec<String>> = BTreeMap::new();
    let mut combined_header = Vec::new();
    let mut sample_names: Vec<String> = Vec::new();

    for (file_idx, input) in inputs.iter().enumerate() {
        let file = File::open(input)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line.map_err(RsomicsError::Io)?;
            if line.starts_with("##") {
                if file_idx == 0 {
                    combined_header.push(line);
                }
                continue;
            }
            if line.starts_with('#') {
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() > 9 {
                    for s in &fields[9..] {
                        sample_names.push(format!("{}:{s}", input.display()));
                    }
                }
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 5 {
                continue;
            }
            let key = (
                fields[0].to_string(),
                fields[1].parse::<u64>().unwrap_or(0),
                fields[3].to_string(),
                fields[4].to_string(),
            );
            all_records.entry(key).or_default().push(line);
        }
    }

    for h in &combined_header {
        writeln!(out, "{h}").map_err(RsomicsError::Io)?;
    }

    let mut count: u64 = 0;
    for (_, lines) in &all_records {
        writeln!(out, "{}", lines[0]).map_err(RsomicsError::Io)?;
        count += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
