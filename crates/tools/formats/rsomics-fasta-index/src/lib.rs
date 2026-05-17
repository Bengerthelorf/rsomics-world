#![allow(clippy::cast_possible_truncation)]

use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct FaiEntry {
    pub name: String,
    pub length: u64,
    pub offset: u64,
    pub line_bases: u64,
    pub line_width: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct FaiIndex {
    pub entries: Vec<FaiEntry>,
}

impl std::fmt::Display for FaiIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for e in &self.entries {
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t{}",
                e.name, e.length, e.offset, e.line_bases, e.line_width
            )?;
        }
        Ok(())
    }
}

pub fn build_index(fasta_path: &Path) -> Result<FaiIndex> {
    let file = File::open(fasta_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", fasta_path.display())))?;
    let mut reader = BufReader::new(file);
    let mut entries = Vec::new();
    let mut line_buf = String::new();
    let mut byte_offset: u64 = 0;

    loop {
        line_buf.clear();
        let n = reader.read_line(&mut line_buf).map_err(RsomicsError::Io)?;
        if n == 0 {
            break;
        }

        if !line_buf.starts_with('>') {
            byte_offset += n as u64;
            continue;
        }

        let name = line_buf[1..]
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        byte_offset += n as u64;
        let seq_start = byte_offset;

        let mut length: u64 = 0;
        let mut line_bases: u64 = 0;
        let mut line_width: u64 = 0;
        let mut first_line = true;

        loop {
            line_buf.clear();
            let ln = reader.read_line(&mut line_buf).map_err(RsomicsError::Io)?;
            if ln == 0 || line_buf.starts_with('>') {
                break;
            }
            let raw_len = ln as u64;
            let bases = line_buf.trim_end_matches(&['\n', '\r'][..]).len() as u64;
            length += bases;
            if first_line {
                line_bases = bases;
                line_width = raw_len;
                first_line = false;
            }
            byte_offset += raw_len;
        }

        entries.push(FaiEntry {
            name,
            length,
            offset: seq_start,
            line_bases,
            line_width,
        });

        if line_buf.starts_with('>') {
            let name = line_buf[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            byte_offset += line_buf.len() as u64;
            let seq_start2 = byte_offset;

            let mut length2: u64 = 0;
            let mut lb2: u64 = 0;
            let mut lw2: u64 = 0;
            let mut first2 = true;

            loop {
                line_buf.clear();
                let ln = reader.read_line(&mut line_buf).map_err(RsomicsError::Io)?;
                if ln == 0 || line_buf.starts_with('>') {
                    break;
                }
                let raw_len = ln as u64;
                let bases = line_buf.trim_end_matches(&['\n', '\r'][..]).len() as u64;
                length2 += bases;
                if first2 {
                    lb2 = bases;
                    lw2 = raw_len;
                    first2 = false;
                }
                byte_offset += raw_len;
            }
            entries.push(FaiEntry {
                name,
                length: length2,
                offset: seq_start2,
                line_bases: lb2,
                line_width: lw2,
            });
        }
    }

    Ok(FaiIndex { entries })
}

pub fn write_index(index: &FaiIndex, path: &Path) -> Result<()> {
    let mut file = File::create(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("creating {}: {e}", path.display())))?;
    write!(file, "{index}").map_err(RsomicsError::Io)?;
    Ok(())
}

pub fn fetch_region(fasta_path: &Path, fai_path: &Path, region: &str) -> Result<Vec<u8>> {
    let index = read_index(fai_path)?;
    let (name, start, end) = parse_region(region)?;

    let entry = index
        .entries
        .iter()
        .find(|e| e.name == name)
        .ok_or_else(|| {
            RsomicsError::InvalidInput(format!("sequence '{name}' not found in index"))
        })?;

    let start = start.unwrap_or(0);
    let end = end.unwrap_or(entry.length as usize);

    if start >= end || start as u64 >= entry.length {
        return Ok(Vec::new());
    }
    let end = end.min(entry.length as usize);

    let mut file = File::open(fasta_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", fasta_path.display())))?;

    let start_line = start as u64 / entry.line_bases;
    let start_col = start as u64 % entry.line_bases;
    let seek_pos = entry.offset + start_line * entry.line_width + start_col;

    file.seek(SeekFrom::Start(seek_pos))
        .map_err(RsomicsError::Io)?;

    let mut result = Vec::with_capacity(end - start);
    let mut buf = [0u8; 8192];
    let mut remaining = (end - start) as u64;

    while remaining > 0 {
        let read_len = buf.len().min(remaining as usize + 1024);
        let n = file.read(&mut buf[..read_len]).map_err(RsomicsError::Io)?;
        if n == 0 {
            break;
        }
        for &b in &buf[..n] {
            if remaining == 0 {
                break;
            }
            if b != b'\n' && b != b'\r' {
                result.push(b);
                remaining -= 1;
            }
        }
    }

    Ok(result)
}

pub fn read_index(path: &Path) -> Result<FaiIndex> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let reader = BufReader::new(file);
    let mut entries = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            return Err(RsomicsError::InvalidInput(format!(
                "invalid .fai line (need 5 fields): {line}"
            )));
        }
        entries.push(FaiEntry {
            name: fields[0].to_string(),
            length: fields[1]
                .parse()
                .map_err(|e| RsomicsError::InvalidInput(format!("bad length in .fai: {e}")))?,
            offset: fields[2]
                .parse()
                .map_err(|e| RsomicsError::InvalidInput(format!("bad offset in .fai: {e}")))?,
            line_bases: fields[3]
                .parse()
                .map_err(|e| RsomicsError::InvalidInput(format!("bad line_bases in .fai: {e}")))?,
            line_width: fields[4]
                .parse()
                .map_err(|e| RsomicsError::InvalidInput(format!("bad line_width in .fai: {e}")))?,
        });
    }

    Ok(FaiIndex { entries })
}

fn parse_region(region: &str) -> Result<(String, Option<usize>, Option<usize>)> {
    if let Some(colon) = region.rfind(':') {
        let name = &region[..colon];
        let range = &region[colon + 1..];
        if let Some(dash) = range.find('-') {
            let start: usize = range[..dash]
                .replace(',', "")
                .parse()
                .map_err(|e| RsomicsError::InvalidInput(format!("bad region start: {e}")))?;
            let end: usize = range[dash + 1..]
                .replace(',', "")
                .parse()
                .map_err(|e| RsomicsError::InvalidInput(format!("bad region end: {e}")))?;
            Ok((name.to_string(), Some(start.saturating_sub(1)), Some(end)))
        } else {
            let start: usize = range
                .replace(',', "")
                .parse()
                .map_err(|e| RsomicsError::InvalidInput(format!("bad region pos: {e}")))?;
            Ok((name.to_string(), Some(start.saturating_sub(1)), Some(start)))
        }
    } else {
        Ok((region.to_string(), None, None))
    }
}
