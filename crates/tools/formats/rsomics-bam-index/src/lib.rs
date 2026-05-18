use std::fs::File;
use std::path::Path;

use noodles::bam;
use noodles::bam::bai;
use rsomics_common::{Result, RsomicsError};

pub fn index_bam(bam_path: &Path) -> Result<()> {
    let bai_path = bam_path.with_extension("bam.bai");

    let mut reader = File::open(bam_path)
        .map(bam::io::Reader::new)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bam_path.display())))?;

    let header = reader
        .read_header()
        .map_err(|e| RsomicsError::InvalidInput(format!("reading header: {e}")))?;

    let index = bam::bai::index(&mut reader, &header)
        .map_err(|e| RsomicsError::InvalidInput(format!("building index: {e}")))?;

    bai::fs::write(&bai_path, &index)
        .map_err(|e| RsomicsError::InvalidInput(format!("writing {}: {e}", bai_path.display())))?;

    Ok(())
}
