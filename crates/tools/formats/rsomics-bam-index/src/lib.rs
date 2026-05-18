use std::path::Path;
use std::process::Command;

use rsomics_common::{Result, RsomicsError};

pub fn index_bam(bam_path: &Path) -> Result<()> {
    let status = Command::new("samtools")
        .args(["index", bam_path.to_str().unwrap_or("-")])
        .status()
        .map_err(|e| {
            RsomicsError::InvalidInput(format!(
                "failed to run samtools index: {e}. Is samtools installed?"
            ))
        })?;

    if !status.success() {
        return Err(RsomicsError::UpstreamError(format!(
            "samtools index exited with {}",
            status.code().unwrap_or(-1)
        )));
    }

    Ok(())
}
