use std::fmt;
use std::fs::File;
use std::io;
use std::path::Path;

use noodles::bam;
use noodles::sam;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

#[derive(Debug, Default, Clone, Serialize)]
pub struct FlagstatCounts {
    pub total: [u64; 2],
    pub secondary: [u64; 2],
    pub supplementary: [u64; 2],
    pub duplicates: [u64; 2],
    pub primary_duplicates: [u64; 2],
    pub mapped: [u64; 2],
    pub primary_mapped: [u64; 2],
    pub paired: [u64; 2],
    pub read1: [u64; 2],
    pub read2: [u64; 2],
    pub properly_paired: [u64; 2],
    pub both_mapped: [u64; 2],
    pub singletons: [u64; 2],
    pub mate_diff_chr: [u64; 2],
    pub mate_diff_chr_mapq5: [u64; 2],
}

fn pct(num: u64, den: u64) -> String {
    if den == 0 {
        "N/A".to_owned()
    } else {
        format!("{:.2}%", num as f64 / den as f64 * 100.0)
    }
}

impl fmt::Display for FlagstatCounts {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let primary = [
            self.total[0] - self.secondary[0] - self.supplementary[0],
            self.total[1] - self.secondary[1] - self.supplementary[1],
        ];
        writeln!(f, "{} + {} in total (QC-passed reads + QC-failed reads)", self.total[0], self.total[1])?;
        writeln!(f, "{} + {} primary", primary[0], primary[1])?;
        writeln!(f, "{} + {} secondary", self.secondary[0], self.secondary[1])?;
        writeln!(f, "{} + {} supplementary", self.supplementary[0], self.supplementary[1])?;
        writeln!(f, "{} + {} duplicates", self.duplicates[0], self.duplicates[1])?;
        writeln!(f, "{} + {} primary duplicates", self.primary_duplicates[0], self.primary_duplicates[1])?;
        writeln!(f, "{} + {} mapped ({} : {})", self.mapped[0], self.mapped[1], pct(self.mapped[0], self.total[0]), pct(self.mapped[1], self.total[1]))?;
        writeln!(f, "{} + {} primary mapped ({} : {})", self.primary_mapped[0], self.primary_mapped[1], pct(self.primary_mapped[0], primary[0]), pct(self.primary_mapped[1], primary[1]))?;
        writeln!(f, "{} + {} paired in sequencing", self.paired[0], self.paired[1])?;
        writeln!(f, "{} + {} read1", self.read1[0], self.read1[1])?;
        writeln!(f, "{} + {} read2", self.read2[0], self.read2[1])?;
        writeln!(f, "{} + {} properly paired ({} : {})", self.properly_paired[0], self.properly_paired[1], pct(self.properly_paired[0], self.paired[0]), pct(self.properly_paired[1], self.paired[1]))?;
        writeln!(f, "{} + {} with itself and mate mapped", self.both_mapped[0], self.both_mapped[1])?;
        writeln!(f, "{} + {} singletons ({} : {})", self.singletons[0], self.singletons[1], pct(self.singletons[0], self.paired[0]), pct(self.singletons[1], self.paired[1]))?;
        writeln!(f, "{} + {} with mate mapped to a different chr", self.mate_diff_chr[0], self.mate_diff_chr[1])?;
        write!(f, "{} + {} with mate mapped to a different chr (mapQ>=5)", self.mate_diff_chr_mapq5[0], self.mate_diff_chr_mapq5[1])?;
        Ok(())
    }
}

pub fn count_bam(path: &Path) -> Result<FlagstatCounts> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header()
        .map_err(|e| RsomicsError::InvalidInput(format!("reading header from {}: {e}", path.display())))?;

    let mut c = FlagstatCounts::default();

    for result in reader.records() {
        let record = result
            .map_err(|e| RsomicsError::InvalidInput(format!("reading record from {}: {e}", path.display())))?;
        tally_record(&record, &header, &mut c)?;
    }

    Ok(c)
}

pub fn count_sam<R: io::BufRead>(reader: R) -> Result<FlagstatCounts> {
    let mut sam_reader = sam::io::Reader::new(reader);
    let header = sam_reader.read_header()
        .map_err(|e| RsomicsError::InvalidInput(format!("reading SAM header: {e}")))?;

    let mut c = FlagstatCounts::default();

    for result in sam_reader.records() {
        let record = result
            .map_err(|e| RsomicsError::InvalidInput(format!("reading SAM record: {e}")))?;
        tally_record(&record, &header, &mut c)?;
    }

    Ok(c)
}

fn tally_record<R: sam::alignment::Record>(record: &R, header: &sam::Header, c: &mut FlagstatCounts) -> Result<()> {
    let flags = record.flags()
        .map_err(|e| RsomicsError::InvalidInput(format!("reading flags: {e}")))?;

    let i = usize::from(flags.is_qc_fail());
    let is_secondary = flags.is_secondary();
    let is_supplementary = flags.is_supplementary();
    let is_primary = !is_secondary && !is_supplementary;

    c.total[i] += 1;

    if is_secondary {
        c.secondary[i] += 1;
    }
    if is_supplementary {
        c.supplementary[i] += 1;
    }
    if flags.is_duplicate() {
        c.duplicates[i] += 1;
        if is_primary {
            c.primary_duplicates[i] += 1;
        }
    }

    let is_mapped = !flags.is_unmapped();
    if is_mapped {
        c.mapped[i] += 1;
        if is_primary {
            c.primary_mapped[i] += 1;
        }
    }

    if flags.is_segmented() {
        c.paired[i] += 1;
        if flags.is_first_segment() {
            c.read1[i] += 1;
        }
        if flags.is_last_segment() {
            c.read2[i] += 1;
        }
        if flags.is_properly_segmented() {
            c.properly_paired[i] += 1;
        }

        let mate_mapped = !flags.is_mate_unmapped();
        if is_mapped && mate_mapped {
            c.both_mapped[i] += 1;

            if is_primary {
                let tid = record.reference_sequence_id(header)
                    .transpose()
                    .map_err(|e| RsomicsError::InvalidInput(format!("tid: {e}")))?;
                let mtid = record.mate_reference_sequence_id(header)
                    .transpose()
                    .map_err(|e| RsomicsError::InvalidInput(format!("mtid: {e}")))?;

                if tid != mtid {
                    c.mate_diff_chr[i] += 1;
                    if let Some(Ok(mapq)) = record.mapping_quality() {
                        if mapq.get() >= 5 {
                            c.mate_diff_chr_mapq5[i] += 1;
                        }
                    }
                }
            }
        }

        if is_mapped && !mate_mapped && is_primary {
            c.singletons[i] += 1;
        }
    }

    Ok(())
}
