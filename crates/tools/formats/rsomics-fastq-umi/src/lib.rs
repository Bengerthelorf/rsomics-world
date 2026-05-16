pub(crate) mod parallel_gz;

use std::path::Path;

use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::{OwnedRecord, open_fastq};
use serde::Serialize;

use crate::parallel_gz::ChunkedWriter;

const CHUNK_RECORDS: usize = 8192;

/// Which read carries the inline 5' UMI (fastp `UMI_LOC_READ1` / `UMI_LOC_READ2`).
///
/// Index-based locations (`index1`/`index2`/`per_index`/`per_read`) need
/// separate index-FASTQ inputs or dual-UMI merge; they are intentionally not
/// implemented until a consumer needs them (the inline read1/read2 case is the
/// dominant UMI workflow and what `umi_tools extract` defaults to).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmiLoc {
    Read1,
    Read2,
}

/// Exact fastp UMI semantics — sourced from fastp `src/umiprocessor.cpp`
/// (`UmiProcessor::process` / `addUmiToName`) and `src/options.h` `UMIOptions`
/// (fastp MIT; reading + citing permitted). `umi_tools` (Smith et al. 2017,
/// doi:10.1101/gr.209601.116, MIT) is the secondary behavioural reference.
///
/// fastp: `umiLength = min(read.length, umi_len)`; the source read is then
/// `trimFront(umiLength + umi_skip)` (UMI + skip removed from the 5' of both
/// seq and qual); the UMI is appended to the read name as
/// `delimiter + [prefix + "_"] + umi`, inserted before the first space if the
/// name has one, else appended. The non-source mate is not trimmed but its
/// name is stamped with the same UMI so the pair stays consistent.
#[derive(Debug, Clone)]
pub struct UmiConfig {
    pub loc: UmiLoc,
    /// fastp `--umi_len`. Bounded per-read by the read length.
    pub len: usize,
    /// fastp `--umi_skip`: extra 5' bases removed after the UMI.
    pub skip: usize,
    /// fastp `--umi_prefix`: joined to the UMI with `_` when non-empty.
    pub prefix: Vec<u8>,
    /// Read-name / UMI separator. fastp default `:`.
    pub delimiter: u8,
}

impl UmiConfig {
    /// Build the name tag fastp would append: `delimiter [prefix _] umi`.
    fn tag(&self, umi: &[u8]) -> Vec<u8> {
        let mut t = Vec::with_capacity(1 + self.prefix.len() + 1 + umi.len());
        t.push(self.delimiter);
        if !self.prefix.is_empty() {
            t.extend_from_slice(&self.prefix);
            t.push(b'_');
        }
        t.extend_from_slice(umi);
        t
    }
}

/// Insert `tag` before the first space in `id`, else append (fastp
/// `addUmiToName` placement).
fn stamp(id: &[u8], tag: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(id.len() + tag.len());
    if let Some(sp) = id.iter().position(|&b| b == b' ') {
        out.extend_from_slice(&id[..sp]);
        out.extend_from_slice(tag);
        out.extend_from_slice(&id[sp..]);
    } else {
        out.extend_from_slice(id);
        out.extend_from_slice(tag);
    }
    out
}

/// Extract the 5' UMI from `src`, trim it off, and stamp the source name.
/// Returns the tag stamped (so PE callers can stamp the mate identically).
///
/// Byte-faithful to fastp 0.20.1: `umi = seq[..min(umi_len, read_len)]`;
/// then fastp `Read::trimFront(len)` does `len = min(length()-1, len)` — it
/// keeps at least one base, so the 5' trim is `min(umi_len + skip,
/// read_len - 1)`, NOT `read_len`. fastp also guards `if(!umi.empty())
/// addUmiToName(...)` and its `trimFront` on a zero-length read throws; a
/// zero-length UMI source read therefore has no defined fastp output, so we
/// fail loud rather than fabricate a record fastp never emits.
///
/// # Errors
///
/// `InvalidInput` if the UMI source read has zero length.
fn apply_umi(src: &mut OwnedRecord, cfg: &UmiConfig) -> Result<Vec<u8>> {
    let read_len = src.seq.len();
    let umi_len = cfg.len.min(read_len);
    if umi_len == 0 {
        return Err(RsomicsError::InvalidInput(
            "UMI source read has zero length; cannot extract a UMI".into(),
        ));
    }
    let umi = src.seq[..umi_len].to_vec();
    let trim = (umi_len + cfg.skip).min(read_len - 1);
    src.seq.drain(..trim);
    src.qual.drain(..trim);
    let tag = cfg.tag(&umi);
    src.id = stamp(&src.id, &tag);
    Ok(tag)
}

#[derive(Debug, Default, Clone, Serialize)]
pub struct UmiReport {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mode: Option<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub input_r1: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub input_r2: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output_r1: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output_r2: Option<String>,
    pub reads_in: u64,
    pub reads_out: u64,
    pub bases_in: u64,
    pub bases_out: u64,
}

struct OwnedPair {
    r1: OwnedRecord,
    r2: OwnedRecord,
}

pub struct Pipeline<'cfg> {
    pub cfg: &'cfg UmiConfig,
    pub compression: i32,
}

impl<'cfg> Pipeline<'cfg> {
    #[must_use]
    pub fn new(cfg: &'cfg UmiConfig, compression: i32) -> Self {
        Self { cfg, compression }
    }

    /// SE: the single read carries the UMI (only `--umi_loc read1` is valid;
    /// `read2` in SE has no second read and is rejected by the CLI).
    ///
    /// # Errors
    ///
    /// Propagates input parse / output write errors.
    pub fn run_se(&self, input: &Path, output: &Path) -> Result<UmiReport> {
        let mut reader = open_fastq(input)?;
        let mut writer = ChunkedWriter::create(output, self.compression)?;
        let mut report = UmiReport {
            mode: Some("SE"),
            input_r1: Some(input.display().to_string()),
            output_r1: Some(output.display().to_string()),
            ..UmiReport::default()
        };
        let mut chunk: Vec<OwnedRecord> = Vec::with_capacity(CHUNK_RECORDS);
        loop {
            chunk.clear();
            while chunk.len() < CHUNK_RECORDS {
                let Some(r) = reader.next() else { break };
                chunk.push(r?);
            }
            if chunk.is_empty() {
                break;
            }
            let out: Vec<(OwnedRecord, u64, u64)> = chunk
                .par_drain(..)
                .map(|mut rec| {
                    let bases_in = rec.seq.len() as u64;
                    apply_umi(&mut rec, self.cfg)?;
                    let bases_out = rec.seq.len() as u64;
                    Ok((rec, bases_in, bases_out))
                })
                .collect::<Result<Vec<_>>>()?;
            for (rec, bin, bout) in out {
                report.reads_in += 1;
                report.reads_out += 1;
                report.bases_in += bin;
                report.bases_out += bout;
                writer.write_record(&rec.id, &rec.seq, &rec.qual)?;
            }
        }
        writer.finalize()?;
        Ok(report)
    }

    /// PE: the UMI comes from the read selected by `cfg.loc`; only that read is
    /// trimmed, both mate names are stamped with the same UMI.
    ///
    /// # Errors
    ///
    /// Propagates input parse / output write errors; errors if the two inputs
    /// have a differing record count.
    pub fn run_pe(&self, in1: &Path, in2: &Path, out1: &Path, out2: &Path) -> Result<UmiReport> {
        let mut r1 = open_fastq(in1)?;
        let mut r2 = open_fastq(in2)?;
        let mut w1 = ChunkedWriter::create(out1, self.compression)?;
        let mut w2 = ChunkedWriter::create(out2, self.compression)?;
        let mut report = UmiReport {
            mode: Some("PE"),
            input_r1: Some(in1.display().to_string()),
            input_r2: Some(in2.display().to_string()),
            output_r1: Some(out1.display().to_string()),
            output_r2: Some(out2.display().to_string()),
            ..UmiReport::default()
        };
        let mut chunk: Vec<OwnedPair> = Vec::with_capacity(CHUNK_RECORDS);
        let mut done = false;
        while !done {
            chunk.clear();
            while chunk.len() < CHUNK_RECORDS {
                match (r1.next(), r2.next()) {
                    (Some(a), Some(b)) => chunk.push(OwnedPair { r1: a?, r2: b? }),
                    (None, None) => {
                        done = true;
                        break;
                    }
                    _ => {
                        return Err(RsomicsError::InvalidInput(
                            "PE input record counts diverge".into(),
                        ));
                    }
                }
            }
            if chunk.is_empty() {
                break;
            }
            let out: Vec<(OwnedPair, u64, u64)> = chunk
                .par_drain(..)
                .map(|mut p| {
                    let bases_in = (p.r1.seq.len() + p.r2.seq.len()) as u64;
                    let tag = match self.cfg.loc {
                        UmiLoc::Read1 => apply_umi(&mut p.r1, self.cfg)?,
                        UmiLoc::Read2 => apply_umi(&mut p.r2, self.cfg)?,
                    };
                    // Stamp the non-source mate's name with the same UMI tag so
                    // the pair keeps a consistent identifier (fastp behaviour).
                    match self.cfg.loc {
                        UmiLoc::Read1 => p.r2.id = stamp(&p.r2.id, &tag),
                        UmiLoc::Read2 => p.r1.id = stamp(&p.r1.id, &tag),
                    }
                    let bases_out = (p.r1.seq.len() + p.r2.seq.len()) as u64;
                    Ok((p, bases_in, bases_out))
                })
                .collect::<Result<Vec<_>>>()?;
            for (p, bin, bout) in out {
                report.reads_in += 2;
                report.reads_out += 2;
                report.bases_in += bin;
                report.bases_out += bout;
                w1.write_record(&p.r1.id, &p.r1.seq, &p.r1.qual)?;
                w2.write_record(&p.r2.id, &p.r2.seq, &p.r2.qual)?;
            }
        }
        w1.finalize()?;
        w2.finalize()?;
        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg(loc: UmiLoc, len: usize, skip: usize, prefix: &str) -> UmiConfig {
        UmiConfig {
            loc,
            len,
            skip,
            prefix: prefix.as_bytes().to_vec(),
            delimiter: b':',
        }
    }

    fn rec(id: &str, seq: &str, qual: &str) -> OwnedRecord {
        OwnedRecord {
            id: id.as_bytes().to_vec(),
            seq: seq.as_bytes().to_vec(),
            qual: qual.as_bytes().to_vec(),
        }
    }

    #[test]
    fn extract_trims_and_stamps_no_space() {
        let mut r = rec("read1", "ACGTACGTAA", "IIIIIIIIII");
        apply_umi(&mut r, &cfg(UmiLoc::Read1, 4, 0, "")).unwrap();
        assert_eq!(r.id, b"read1:ACGT");
        assert_eq!(r.seq, b"ACGTAA");
        assert_eq!(r.qual, b"IIIIII");
    }

    #[test]
    fn stamp_inserts_before_first_space() {
        let mut r = rec("read1 1:N:0", "TTTTGGGG", "IIIIIIII");
        apply_umi(&mut r, &cfg(UmiLoc::Read1, 4, 0, "")).unwrap();
        assert_eq!(r.id, b"read1:TTTT 1:N:0");
        assert_eq!(r.seq, b"GGGG");
    }

    #[test]
    fn skip_removes_extra_bases_after_umi() {
        let mut r = rec("r", "AACCGGTT", "IIIIIIII");
        apply_umi(&mut r, &cfg(UmiLoc::Read1, 2, 2, "")).unwrap();
        assert_eq!(r.id, b"r:AA");
        assert_eq!(r.seq, b"GGTT");
        assert_eq!(r.qual, b"IIII");
    }

    #[test]
    fn prefix_joined_with_underscore() {
        let mut r = rec("r", "ACGTACGT", "IIIIIIII");
        apply_umi(&mut r, &cfg(UmiLoc::Read1, 4, 0, "UMI")).unwrap();
        assert_eq!(r.id, b"r:UMI_ACGT");
    }

    /// fastp `Read::trimFront` clamps to `length()-1`, so a read shorter than
    /// the requested trim keeps its last base — it is never emptied.
    #[test]
    fn umi_len_clamped_keeps_last_base() {
        let mut r = rec("r", "ACG", "III");
        apply_umi(&mut r, &cfg(UmiLoc::Read1, 8, 0, "")).unwrap();
        // umi = min(8,3)=3 → "ACG"; trim = min(3+0, 3-1) = 2 → "G"/"I" remain.
        assert_eq!(r.id, b"r:ACG");
        assert_eq!(r.seq, b"G");
        assert_eq!(r.qual, b"I");
    }

    #[test]
    fn skip_overrun_keeps_last_base() {
        let mut r = rec("r", "AACCGGTT", "IIIIIIIJ");
        apply_umi(&mut r, &cfg(UmiLoc::Read1, 4, 20, "")).unwrap();
        // umi "AACC"; trim = min(4+20, 8-1) = 7 → last base kept.
        assert_eq!(r.id, b"r:AACC");
        assert_eq!(r.seq, b"T");
        assert_eq!(r.qual, b"J");
    }

    #[test]
    fn exact_consume_keeps_last_base() {
        let mut r = rec("r", "AACCGGTT", "IIIIIIIJ");
        apply_umi(&mut r, &cfg(UmiLoc::Read1, 4, 4, "")).unwrap();
        // umi "AACC"; trim = min(4+4, 8-1) = 7 → last base kept.
        assert_eq!(r.seq, b"T");
        assert_eq!(r.qual, b"J");
    }

    #[test]
    fn one_base_read_kept_and_stamped() {
        let mut r = rec("r", "A", "I");
        apply_umi(&mut r, &cfg(UmiLoc::Read1, 1, 0, "")).unwrap();
        // umi "A"; trim = min(1+0, 1-1) = 0 → the 1 base is kept.
        assert_eq!(r.id, b"r:A");
        assert_eq!(r.seq, b"A");
        assert_eq!(r.qual, b"I");
    }

    #[test]
    fn empty_source_read_errors() {
        let mut r = rec("r", "", "");
        assert!(apply_umi(&mut r, &cfg(UmiLoc::Read1, 4, 0, "")).is_err());
    }
}
