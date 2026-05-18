use std::cmp::Ordering;
use std::fs::File;
use std::path::Path;

use noodles::bam;
use rsomics_common::{Result, RsomicsError};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    Coordinate,
    Name,
}

fn tid(r: &bam::Record) -> Option<usize> {
    r.reference_sequence_id().transpose().ok().flatten()
}

fn pos(r: &bam::Record) -> Option<usize> {
    r.alignment_start()
        .transpose()
        .ok()
        .flatten()
        .map(|p| p.get())
}

pub fn sort_bam(input: &Path, output: &Path, order: SortOrder) -> Result<()> {
    let file = File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader
        .read_header()
        .map_err(|e| RsomicsError::InvalidInput(format!("reading header: {e}")))?;

    let mut records: Vec<bam::Record> = Vec::new();
    for result in reader.records() {
        let record =
            result.map_err(|e| RsomicsError::InvalidInput(format!("reading record: {e}")))?;
        records.push(record);
    }

    match order {
        SortOrder::Coordinate => {
            records.sort_by(|a, b| match (tid(a), tid(b)) {
                (None, None) => Ordering::Equal,
                (None, Some(_)) => Ordering::Greater,
                (Some(_), None) => Ordering::Less,
                (Some(ta), Some(tb)) => ta.cmp(&tb).then_with(|| pos(a).cmp(&pos(b))),
            });
        }
        SortOrder::Name => {
            records.sort_by(|a, b| a.name().cmp(&b.name()));
        }
    }

    let out_file = File::create(output)
        .map_err(|e| RsomicsError::InvalidInput(format!("creating {}: {e}", output.display())))?;
    let mut writer = bam::io::Writer::new(out_file);
    writer
        .write_header(&header)
        .map_err(|e| RsomicsError::InvalidInput(format!("writing header: {e}")))?;

    for record in &records {
        writer
            .write_record(&header, record)
            .map_err(|e| RsomicsError::InvalidInput(format!("writing record: {e}")))?;
    }

    Ok(())
}
