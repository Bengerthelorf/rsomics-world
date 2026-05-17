use std::fs::File;
use std::path::Path;

use noodles::bam;
use noodles::bam::bai;
use noodles::sam::header::record::value::map::ReferenceSequence;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct RefStats {
    pub name: String,
    pub length: usize,
    pub mapped: u64,
    pub unmapped: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct IdxStats {
    pub refs: Vec<RefStats>,
    pub unmapped_no_ref: u64,
}

impl std::fmt::Display for IdxStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for r in &self.refs {
            writeln!(f, "{}\t{}\t{}\t{}", r.name, r.length, r.mapped, r.unmapped)?;
        }
        write!(f, "*\t0\t0\t{}", self.unmapped_no_ref)?;
        Ok(())
    }
}

pub fn idxstats(bam_path: &Path) -> Result<IdxStats> {
    let bai_path = bam_path.with_extension("bam.bai");
    let alt_bai = bam_path.with_file_name(format!(
        "{}.bai",
        bam_path.file_name().unwrap_or_default().to_string_lossy()
    ));

    let index_path = if bai_path.exists() {
        bai_path
    } else if alt_bai.exists() {
        alt_bai
    } else {
        return Err(RsomicsError::InvalidInput(format!(
            "no BAI index found for {} (tried {} and {})",
            bam_path.display(),
            bai_path.display(),
            alt_bai.display()
        )));
    };

    let index = bai::read(&index_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("reading index {}: {e}", index_path.display())))?;

    let file = File::open(bam_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bam_path.display())))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader
        .read_header()
        .map_err(|e| RsomicsError::InvalidInput(format!("reading header: {e}")))?;

    let ref_seqs = header.reference_sequences();
    let mut stats = IdxStats {
        refs: Vec::with_capacity(ref_seqs.len()),
        unmapped_no_ref: 0,
    };

    for (i, (name, map)) in ref_seqs.iter().enumerate() {
        let length = map.length();
        let (mapped, unmapped) = if let Some(ref_idx) = index.indices().get(i) {
            let m = ref_idx
                .metadata()
                .map_or(0, |md| md.mapped_record_count());
            let u = ref_idx
                .metadata()
                .map_or(0, |md| md.unmapped_record_count());
            (m, u)
        } else {
            (0, 0)
        };
        stats.refs.push(RefStats {
            name: name.to_string(),
            length,
            mapped,
            unmapped,
        });
    }

    if let Some(n_no_coor) = index.unplaced_unmapped_record_count() {
        stats.unmapped_no_ref = n_no_coor;
    }

    Ok(stats)
}
