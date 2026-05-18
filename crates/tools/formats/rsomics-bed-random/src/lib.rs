use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rsomics_common::{Result, RsomicsError};

pub fn random_bed(
    genome_path: &Path,
    n: u64,
    length: u64,
    seed: u64,
    output: &mut dyn Write,
) -> Result<()> {
    let genome = load_genome(genome_path)?;
    let chroms: Vec<(&String, &u64)> = genome.iter().collect();
    if chroms.is_empty() {
        return Err(RsomicsError::InvalidInput("empty genome file".into()));
    }
    let total_len: u64 = chroms.iter().map(|(_, &l)| l).sum();
    let mut rng = StdRng::seed_from_u64(seed);
    let mut out = BufWriter::with_capacity(64 * 1024, output);

    for _ in 0..n {
        let mut pos_in_genome = rng.gen_range(0..total_len.saturating_sub(length));
        let mut chrom_name = "";
        for (name, &clen) in &chroms {
            if pos_in_genome < clen {
                chrom_name = name;
                break;
            }
            pos_in_genome -= clen;
        }
        let start = pos_in_genome;
        let end = start + length;
        writeln!(out, "{chrom_name}\t{start}\t{end}").map_err(RsomicsError::Io)?;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

fn load_genome(path: &Path) -> Result<HashMap<String, u64>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();
    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 2 {
            let len: u64 = fields[1]
                .parse()
                .map_err(|e| RsomicsError::InvalidInput(format!("bad length: {e}")))?;
            map.insert(fields[0].to_string(), len);
        }
    }
    Ok(map)
}
