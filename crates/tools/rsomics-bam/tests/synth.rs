//! Synthetic BAM construction used by the integration test suite. A
//! tiny BAM is generated at test time rather than checked into git so
//! the fixture is reproducible from source.

use std::path::Path;

use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::record::Record;
use rust_htslib::bam::{Format, Writer};

/// Write `n` unaligned (flag = 4) reads to `path`. Each record carries
/// a fixed 10-base sequence and a fixed 10-base quality. Returns the
/// number of records actually written (== `n` on success).
///
/// # Panics
///
/// Panics if the BAM writer cannot be opened at `path` or a record fails
/// to serialise — tests rely on this hard-fail since either failure
/// indicates a broken fixture setup, not a condition to recover from.
#[must_use]
pub fn write_unaligned_bam(path: &Path, n: usize) -> usize {
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "unsorted");
    header.push_record(&hd);

    let mut writer = Writer::from_path(path, &header, Format::Bam).expect("create BAM writer");
    let seq = b"ACGTACGTAC";
    let qual = [30u8; 10];
    for i in 0..n {
        let mut rec = Record::new();
        let qname = format!("read_{i:06}");
        rec.set(qname.as_bytes(), None, seq, &qual);
        rec.set_unmapped();
        writer.write(&rec).expect("write record");
    }
    drop(writer);
    n
}
