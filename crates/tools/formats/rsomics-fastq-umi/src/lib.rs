use std::path::Path;

use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::{OwnedRecord, open_fastq};
use serde::Serialize;

use rsomics_fqgz::ChunkedWriter;

const CHUNK_RECORDS: usize = 8192;

// fastp 0.20.1 --umi_loc set: Read1/Read2 trim 5' seq; Index1/Index2 use read-name index field (no trim); PerIndex = firstIndex_lastIndex; PerRead = umi1_umi2, both trimmed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmiLoc {
    Read1,
    Read2,
    Index1,
    Index2,
    PerIndex,
    PerRead,
}

// fastp 0.20.1 Read::lastIndex (src/read.cpp): backward from len-3, returns everything after last ':'/'+''; empty when name < 5 bytes or no delimiter.
fn last_index(name: &[u8]) -> Vec<u8> {
    let len = name.len();
    if len < 5 {
        return Vec::new();
    }
    for i in (0..=len - 3).rev() {
        if name[i] == b':' || name[i] == b'+' {
            return name[i + 1..len].to_vec();
        }
    }
    Vec::new()
}

// fastp 0.20.1 Read::firstIndex (src/read.cpp): backward from len-3; '+' sets field end to index-1 (dual-index split), last ':' ends scan; returns the ':'..'+' (or ':'..end) field.
fn first_index(name: &[u8]) -> Vec<u8> {
    let len = name.len();
    if len < 5 {
        return Vec::new();
    }
    let mut end = len;
    for i in (0..=len - 3).rev() {
        if name[i] == b'+' {
            end = i.saturating_sub(1);
        }
        if name[i] == b':' {
            // fastp substr(i+1, end-i) == name[i+1 ..= end] → half-open end+1.
            let stop = (end + 1).min(len);
            let start = (i + 1).min(stop);
            return name[start..stop].to_vec();
        }
    }
    Vec::new()
}

// fastp src/umiprocessor.cpp (MIT): umiLength = min(read.length, umi_len); source read trimFront(umiLength+skip); tag = delimiter[+prefix+"_"]+umi inserted before first space or appended; non-source mate stamped but not trimmed.
#[derive(Debug, Clone)]
pub struct UmiConfig {
    pub loc: UmiLoc,
    pub len: usize,      // fastp --umi_len
    pub skip: usize,     // fastp --umi_skip: extra 5' bases removed after the UMI
    pub prefix: Vec<u8>, // fastp --umi_prefix: joined to the UMI with _ when non-empty
    pub delimiter: u8,   // read-name / UMI separator; fastp default :
}

impl UmiConfig {
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

// fastp 0.20.1 trimFront clamps to length()-1 (keeps ≥1 base); trim = min(umi_len+skip, read_len-1). Zero-length read has no defined fastp output — fail loud.
fn take_seq_umi(src: &mut OwnedRecord, cfg: &UmiConfig) -> Result<Vec<u8>> {
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
    Ok(umi)
}

// fastp 0.20.1 UmiProcessor::process: if(!umi.empty()) guard means a missing index field is a pass-through, not an error.
fn process(
    rec: &mut OwnedRecord,
    mut mate: Option<&mut OwnedRecord>,
    cfg: &UmiConfig,
) -> Result<()> {
    let umi: Vec<u8> = match cfg.loc {
        UmiLoc::Read1 => take_seq_umi(rec, cfg)?,
        UmiLoc::Read2 => {
            let m = mate.as_deref_mut().ok_or_else(|| {
                RsomicsError::ConfigError("--umi_loc read2 requires PE input".into())
            })?;
            take_seq_umi(m, cfg)?
        }
        UmiLoc::Index1 => first_index(&rec.id),
        UmiLoc::Index2 => {
            let m = mate.as_deref().ok_or_else(|| {
                RsomicsError::ConfigError("--umi_loc index2 requires PE input".into())
            })?;
            last_index(&m.id)
        }
        UmiLoc::PerIndex => {
            let mut u = first_index(&rec.id);
            if let Some(m) = mate.as_deref() {
                u.push(b'_');
                u.extend_from_slice(&last_index(&m.id));
            }
            u
        }
        UmiLoc::PerRead => {
            let mut u = take_seq_umi(rec, cfg)?;
            if let Some(m) = mate.as_deref_mut() {
                u.push(b'_');
                u.extend_from_slice(&take_seq_umi(m, cfg)?);
            }
            u
        }
    };
    if umi.is_empty() {
        return Ok(());
    }
    let tag = cfg.tag(&umi);
    rec.id = stamp(&rec.id, &tag);
    if let Some(m) = mate {
        m.id = stamp(&m.id, &tag);
    }
    Ok(())
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
                    process(&mut rec, None, self.cfg)?;
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
                    process(&mut p.r1, Some(&mut p.r2), self.cfg)?;
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
        process(&mut r, None, &cfg(UmiLoc::Read1, 4, 0, "")).unwrap();
        assert_eq!(r.id, b"read1:ACGT");
        assert_eq!(r.seq, b"ACGTAA");
        assert_eq!(r.qual, b"IIIIII");
    }

    #[test]
    fn stamp_inserts_before_first_space() {
        let mut r = rec("read1 1:N:0", "TTTTGGGG", "IIIIIIII");
        process(&mut r, None, &cfg(UmiLoc::Read1, 4, 0, "")).unwrap();
        assert_eq!(r.id, b"read1:TTTT 1:N:0");
        assert_eq!(r.seq, b"GGGG");
    }

    #[test]
    fn skip_removes_extra_bases_after_umi() {
        let mut r = rec("r", "AACCGGTT", "IIIIIIII");
        process(&mut r, None, &cfg(UmiLoc::Read1, 2, 2, "")).unwrap();
        assert_eq!(r.id, b"r:AA");
        assert_eq!(r.seq, b"GGTT");
        assert_eq!(r.qual, b"IIII");
    }

    #[test]
    fn prefix_joined_with_underscore() {
        let mut r = rec("r", "ACGTACGT", "IIIIIIII");
        process(&mut r, None, &cfg(UmiLoc::Read1, 4, 0, "UMI")).unwrap();
        assert_eq!(r.id, b"r:UMI_ACGT");
    }

    // fastp Read::trimFront clamps to length()-1 — a short read always keeps its last base.
    #[test]
    fn umi_len_clamped_keeps_last_base() {
        let mut r = rec("r", "ACG", "III");
        process(&mut r, None, &cfg(UmiLoc::Read1, 8, 0, "")).unwrap();
        assert_eq!(r.id, b"r:ACG");
        assert_eq!(r.seq, b"G");
        assert_eq!(r.qual, b"I");
    }

    #[test]
    fn skip_overrun_keeps_last_base() {
        let mut r = rec("r", "AACCGGTT", "IIIIIIIJ");
        process(&mut r, None, &cfg(UmiLoc::Read1, 4, 20, "")).unwrap();
        assert_eq!(r.id, b"r:AACC");
        assert_eq!(r.seq, b"T");
        assert_eq!(r.qual, b"J");
    }

    #[test]
    fn exact_consume_keeps_last_base() {
        let mut r = rec("r", "AACCGGTT", "IIIIIIIJ");
        process(&mut r, None, &cfg(UmiLoc::Read1, 4, 4, "")).unwrap();
        assert_eq!(r.seq, b"T");
        assert_eq!(r.qual, b"J");
    }

    #[test]
    fn one_base_read_kept_and_stamped() {
        let mut r = rec("r", "A", "I");
        process(&mut r, None, &cfg(UmiLoc::Read1, 1, 0, "")).unwrap();
        assert_eq!(r.id, b"r:A");
        assert_eq!(r.seq, b"A");
        assert_eq!(r.qual, b"I");
    }

    #[test]
    fn empty_source_read_errors() {
        let mut r = rec("r", "", "");
        assert!(process(&mut r, None, &cfg(UmiLoc::Read1, 4, 0, "")).is_err());
    }

    #[test]
    fn index_parsing_matches_fastp_readcpp() {
        // single index: first == last == the trailing field
        assert_eq!(first_index(b"R1 1:N:0:ATCACG"), b"ATCACG");
        assert_eq!(last_index(b"R1 1:N:0:ATCACG"), b"ATCACG");
        // dual index: first = before '+', last = after '+'
        assert_eq!(first_index(b"R1 1:N:0:ATCACG+TGGTCA"), b"ATCACG");
        assert_eq!(last_index(b"R1 1:N:0:ATCACG+TGGTCA"), b"TGGTCA");
        // no delimiter / too short → empty
        assert_eq!(first_index(b"abcd"), b"");
        assert_eq!(last_index(b"readname"), b"");
    }

    #[test]
    fn index1_stamps_from_header_no_trim() {
        let mut r = rec("read 1:N:0:ACGTAA", "TTTTGGGG", "IIIIIIII");
        process(&mut r, None, &cfg(UmiLoc::Index1, 0, 0, "")).unwrap();
        assert_eq!(r.id, b"read:ACGTAA 1:N:0:ACGTAA");
        assert_eq!(r.seq, b"TTTTGGGG"); // index mode never trims
        assert_eq!(r.qual, b"IIIIIIII");
    }

    #[test]
    fn index_missing_field_is_passthrough_not_error() {
        let mut r = rec("plainname", "ACGT", "IIII");
        process(&mut r, None, &cfg(UmiLoc::Index1, 0, 0, "")).unwrap();
        assert_eq!(r.id, b"plainname"); // empty UMI ⇒ fastp !empty guard ⇒ no stamp
        assert_eq!(r.seq, b"ACGT");
    }

    #[test]
    fn per_index_pe_merges_first_and_last() {
        let mut r1 = rec("p 1:N:0:AAA+CCC", "GGGG", "IIII");
        let mut r2 = rec("p 2:N:0:AAA+CCC", "TTTT", "FFFF");
        process(&mut r1, Some(&mut r2), &cfg(UmiLoc::PerIndex, 0, 0, "")).unwrap();
        assert_eq!(r1.id, b"p:AAA_CCC 1:N:0:AAA+CCC");
        assert_eq!(r2.id, b"p:AAA_CCC 2:N:0:AAA+CCC");
        assert_eq!(r1.seq, b"GGGG");
        assert_eq!(r2.seq, b"TTTT");
    }

    #[test]
    fn per_read_pe_merges_both_seq_umis_and_trims_both() {
        let mut r1 = rec("p 1", "AACCGGGG", "IIIIIIII");
        let mut r2 = rec("p 2", "TTGGCCCC", "FFFFFFFF");
        process(&mut r1, Some(&mut r2), &cfg(UmiLoc::PerRead, 4, 0, "")).unwrap();
        assert_eq!(r1.id, b"p:AACC_TTGG 1");
        assert_eq!(r2.id, b"p:AACC_TTGG 2");
        assert_eq!(r1.seq, b"GGGG");
        assert_eq!(r2.seq, b"CCCC");
    }

    #[test]
    fn read2_and_index2_require_pe() {
        let mut r = rec("r 1:N:0:ACGT", "ACGTACGT", "IIIIIIII");
        assert!(process(&mut r, None, &cfg(UmiLoc::Read2, 4, 0, "")).is_err());
        let mut r2 = rec("r 1:N:0:ACGT", "ACGTACGT", "IIIIIIII");
        assert!(process(&mut r2, None, &cfg(UmiLoc::Index2, 0, 0, "")).is_err());
    }
}
