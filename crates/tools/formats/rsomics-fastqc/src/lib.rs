pub mod modules;
pub mod report;

use std::path::Path;

use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};

use modules::{QcModule, Record, default_modules};

pub fn analyze(path: &Path) -> Result<Vec<Box<dyn QcModule>>> {
    let mut reader = parse_fastx_file(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("opening {}: {e}", path.display())))?;

    let display = path.file_name().map_or_else(
        || path.display().to_string(),
        |n| n.to_string_lossy().into(),
    );
    let mut mods = default_modules(&display);

    let mut n_records: u64 = 0;
    while let Some(record) = reader.next() {
        let rec = record
            .map_err(|e| RsomicsError::InvalidInput(format!("parsing {}: {e}", path.display())))?;
        let seq_cow = rec.seq();
        let seq: &[u8] = &seq_cow;
        let qual = rec.qual().ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "{}: record {} has no quality line — not a FASTQ",
                path.display(),
                n_records + 1
            ))
        })?;
        let id = rec.id();
        let r = Record { id, seq, qual };
        for m in &mut mods {
            m.process(&r);
        }
        n_records += 1;
    }

    if n_records == 0 {
        return Err(RsomicsError::InvalidInput(format!(
            "{} contained no FASTQ records",
            path.display()
        )));
    }

    for m in &mut mods {
        m.finalize();
    }
    Ok(mods)
}
