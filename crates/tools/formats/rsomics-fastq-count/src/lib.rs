use std::path::Path;

use rsomics_common::Result;
use rsomics_seqio::open_fastq;

pub fn count(input: &Path) -> Result<u64> {
    let mut reader = open_fastq(input)?;
    let mut n: u64 = 0;
    for result in reader.by_ref() {
        result?;
        n += 1;
    }
    Ok(n)
}
