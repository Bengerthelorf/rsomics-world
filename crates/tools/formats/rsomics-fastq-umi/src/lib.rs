use std::path::Path;

use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::{OwnedRecord, open_fastq};
use serde::Serialize;

use rsomics_fqgz::ChunkedWriter;

const CHUNK_RECORDS: usize = 8192;

/// Where the UMI is taken from — the full fastp 0.20.1 `--umi_loc` set.
///
/// `Read1`/`Read2`: 5' of that read's sequence (the read is trimmed).
/// `Index1`/`Index2`: the read-name trailing index field (`firstIndex` of R1
/// / `lastIndex` of R2), no trim. `PerIndex`: `firstIndex(R1) "_" lastIndex(R2)`.
/// `PerRead`: 5' of both reads' sequences merged `umi1 "_" umi2`, both trimmed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmiLoc {
    Read1,
    Read2,
    Index1,
    Index2,
    PerIndex,
    PerRead,
}

/// fastp 0.20.1 `Read::lastIndex` (src/read.cpp), verbatim: scan the name
/// backward from `len-3` for the last `:`/`+`; return everything after it.
/// Empty when the name is < 5 bytes or has no delimiter. The C `substr`
/// length is clamped to the end, which the slice replicates.
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

/// fastp 0.20.1 `Read::firstIndex` (src/read.cpp), verbatim: backward from
/// `len-3`, a `+` sets the field end to its index-1 (dual-index split), the
/// last `:` ends the scan; return the `:`..`+` (or `:`..end) field. Bounds
/// are clamped exactly as the C `substr` would.
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

/// Take the 5' UMI from `src`'s sequence and trim it (+ skip) off seq & qual.
///
/// Byte-faithful to fastp 0.20.1: `umi = seq[..min(umi_len, read_len)]`; fastp
/// `Read::trimFront` clamps to `length()-1` (keeps ≥1 base) so the trim is
/// `min(umi_len + skip, read_len - 1)`. fastp's `trimFront` on a zero-length
/// read throws (no defined output), so a zero-length seq-UMI source fails loud
/// rather than fabricate a record fastp never emits.
///
/// # Errors
///
/// `InvalidInput` if the seq-UMI source read has zero length.
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

/// Apply the UMI transform to a record (PE: pass `mate`) per `cfg.loc`,
/// mirroring fastp 0.20.1 `UmiProcessor::process`: build the UMI for the
/// location, trim only the sequence-UMI source read(s), then — iff the UMI is
/// non-empty (fastp's `if(!umi.empty())` guard, so a missing index field is a
/// pass-through, not an error) — stamp the same tag into both names.
///
/// # Errors
///
/// `InvalidInput` if a sequence-UMI source read has zero length, or the
/// location needs a mate that is absent (rejected earlier by the CLI).
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

    /// fastp `Read::trimFront` clamps to `length()-1`, so a read shorter than
    /// the requested trim keeps its last base — it is never emptied.
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
        assert_eq!(r1.seq, b"GGGG"); // index mode: no trim either mate
        assert_eq!(r2.seq, b"TTTT");
    }

    #[test]
    fn per_read_pe_merges_both_seq_umis_and_trims_both() {
        let mut r1 = rec("p 1", "AACCGGGG", "IIIIIIII");
        let mut r2 = rec("p 2", "TTGGCCCC", "FFFFFFFF");
        process(&mut r1, Some(&mut r2), &cfg(UmiLoc::PerRead, 4, 0, "")).unwrap();
        assert_eq!(r1.id, b"p:AACC_TTGG 1");
        assert_eq!(r2.id, b"p:AACC_TTGG 2");
        assert_eq!(r1.seq, b"GGGG"); // both trimmed by 4
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
