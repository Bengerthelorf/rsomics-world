use std::io::{self, BufRead, BufReader, Read, Write};

use rsomics_common::{Context, Result, RsomicsError};

use crate::interval::{Interval, Strand};

/// Read BED3/BED6 lines. Empty lines and `#`-prefixed comment lines are
/// skipped, matching bedtools behaviour.
#[allow(clippy::missing_errors_doc)]
pub fn read<R: Read>(r: R) -> Result<Vec<Interval>> {
    let mut out = Vec::new();
    for (lineno, line) in BufReader::new(r).lines().enumerate() {
        let line = line.map_err(RsomicsError::Io)?;
        let trimmed = line.trim_end();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let iv = parse_line(trimmed)
            .map_err(|e| RsomicsError::InvalidInput(format!("BED line {}: {e}", lineno + 1)))?;
        out.push(iv);
    }
    Ok(out)
}

fn parse_line(s: &str) -> std::result::Result<Interval, String> {
    let mut fields = s.split('\t');
    let chrom = fields.next().ok_or("missing chrom")?;
    let start_s = fields.next().ok_or("missing start")?;
    let end_s = fields.next().ok_or("missing end")?;
    let start: u64 = start_s
        .parse()
        .map_err(|e| format!("bad start {start_s:?}: {e}"))?;
    let end: u64 = end_s
        .parse()
        .map_err(|e| format!("bad end {end_s:?}: {e}"))?;
    let _name = fields.next();
    let _score = fields.next();
    let strand = fields
        .next()
        .and_then(|s| s.as_bytes().first().copied())
        .and_then(Strand::from_byte);
    let mut iv = Interval::new(chrom, start, end).map_err(|e| e.to_string())?;
    iv.strand = strand;
    Ok(iv)
}

#[allow(clippy::missing_errors_doc)]
pub fn write_bed3<W: Write, I: IntoIterator<Item = Interval>>(mut w: W, ivs: I) -> Result<()> {
    for iv in ivs {
        writeln!(w, "{}\t{}\t{}", iv.chrom, iv.start, iv.end).rs_context("writing BED3 record")?;
    }
    Ok(())
}

/// Write as BED6. `name`/`score` emit `.`/`0` (bedtools placeholders);
/// strand emits `+`/`-`/`.`.
#[allow(clippy::missing_errors_doc)]
pub fn write_bed6<W: Write, I: IntoIterator<Item = Interval>>(mut w: W, ivs: I) -> Result<()> {
    for iv in ivs {
        let strand_byte = iv.strand.map_or('.', |s| s.as_byte() as char);
        writeln!(
            w,
            "{}\t{}\t{}\t.\t0\t{strand_byte}",
            iv.chrom, iv.start, iv.end
        )
        .rs_context("writing BED6 record")?;
    }
    Ok(())
}

#[allow(clippy::missing_errors_doc)]
pub fn read_bytes(bytes: &[u8]) -> Result<Vec<Interval>> {
    read(io::Cursor::new(bytes))
}

/// Single-pass merge of overlapping/touching intervals over an input that
/// is already sorted by chrom then start — the `bedtools merge` contract.
/// Out-of-order input (a smaller start on the active chrom, or a chrom
/// that already closed reappearing) fails loud rather than emitting a
/// wrong result. `track`/`browser`/`#` preamble lines are skipped like
/// bedtools. chrom is treated as opaque bytes, written through verbatim
/// (bedtools never UTF-8-validates a contig name). No per-line `String`,
/// no full materialisation: one reused line buffer, the running cluster,
/// and a chrom buffer that reallocates only on chrom change.
#[allow(clippy::missing_errors_doc)]
pub fn merge_sorted<R: Read, W: Write>(r: R, mut w: W) -> Result<()> {
    let mut rdr = BufReader::new(r);
    let mut line: Vec<u8> = Vec::with_capacity(256);
    let mut chrom: Vec<u8> = Vec::with_capacity(32);
    // Closed chroms, for the reappearing-chrom sortedness check. One entry
    // per distinct chrom (tens at most) — a linear scan beats hashing here.
    let mut closed: Vec<Vec<u8>> = Vec::new();
    let mut have = false;
    let (mut cstart, mut cend) = (0_u64, 0_u64);
    let mut lineno = 0_usize;

    loop {
        line.clear();
        if rdr.read_until(b'\n', &mut line).map_err(RsomicsError::Io)? == 0 {
            break;
        }
        lineno += 1;
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        if line.is_empty()
            || line[0] == b'#'
            || line.starts_with(b"track")
            || line.starts_with(b"browser")
        {
            continue;
        }
        let (lc, ls, le) = parse_bed3_bytes(&line)
            .map_err(|e| RsomicsError::InvalidInput(format!("BED line {lineno}: {e}")))?;

        if have && lc == chrom.as_slice() {
            if ls < cstart {
                return Err(RsomicsError::InvalidInput(format!(
                    "BED line {lineno}: input not sorted (start {ls} < {cstart} on same chrom) — \
                     sort with `sort -k1,1 -k2,2n` first"
                )));
            }
            if ls <= cend {
                cend = cend.max(le);
                continue;
            }
            emit(&mut w, &chrom, cstart, cend)?;
            cstart = ls;
            cend = le;
        } else {
            if have {
                emit(&mut w, &chrom, cstart, cend)?;
                closed.push(chrom.clone());
            }
            if closed.iter().any(|c| c.as_slice() == lc) {
                return Err(RsomicsError::InvalidInput(format!(
                    "BED line {lineno}: input not sorted (chromosome {} reappears after it \
                     closed) — sort with `sort -k1,1 -k2,2n` first",
                    String::from_utf8_lossy(lc)
                )));
            }
            chrom.clear();
            chrom.extend_from_slice(lc);
            cstart = ls;
            cend = le;
            have = true;
        }
    }
    if have {
        emit(&mut w, &chrom, cstart, cend)?;
    }
    Ok(())
}

fn emit<W: Write>(w: &mut W, chrom: &[u8], start: u64, end: u64) -> Result<()> {
    w.write_all(chrom).rs_context("writing merged BED")?;
    writeln!(w, "\t{start}\t{end}").rs_context("writing merged BED")?;
    Ok(())
}

fn parse_bed3_bytes(s: &[u8]) -> std::result::Result<(&[u8], u64, u64), String> {
    let mut it = s.split(|&c| c == b'\t');
    let chrom = it.next().ok_or("missing chrom")?;
    let start = parse_u64(it.next().ok_or("missing start")?)?;
    let end = parse_u64(it.next().ok_or("missing end")?)?;
    if start >= end {
        return Err(format!("empty or inverted interval: {start} >= {end}"));
    }
    Ok((chrom, start, end))
}

fn parse_u64(b: &[u8]) -> std::result::Result<u64, String> {
    if b.is_empty() {
        return Err("empty integer field".into());
    }
    let mut n: u64 = 0;
    for &c in b {
        let d = c.wrapping_sub(b'0');
        if d > 9 {
            return Err(format!("bad integer {:?}", String::from_utf8_lossy(b)));
        }
        n = n
            .checked_mul(10)
            .and_then(|n| n.checked_add(u64::from(d)))
            .ok_or_else(|| format!("integer overflows u64: {:?}", String::from_utf8_lossy(b)))?;
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_bed3_three_records() {
        let input = "chr1\t100\t200\nchr1\t300\t400\nchr2\t50\t150\n";
        let ivs = read_bytes(input.as_bytes()).unwrap();
        assert_eq!(ivs.len(), 3);
        assert_eq!(ivs[0].chrom, "chr1");
        assert_eq!((ivs[0].start, ivs[0].end), (100, 200));
        assert!(ivs[0].strand.is_none());
    }

    #[test]
    fn read_bed6_captures_strand() {
        let input = "chr1\t100\t200\tgene_a\t0\t+\nchr1\t300\t400\tgene_b\t0\t-\n";
        let ivs = read_bytes(input.as_bytes()).unwrap();
        assert_eq!(ivs[0].strand, Some(Strand::Forward));
        assert_eq!(ivs[1].strand, Some(Strand::Reverse));
    }

    #[test]
    fn skip_comments_and_blank_lines() {
        let input = "# header\n\nchr1\t100\t200\n# trailing comment\nchr1\t300\t400\n";
        let ivs = read_bytes(input.as_bytes()).unwrap();
        assert_eq!(ivs.len(), 2);
    }

    #[test]
    fn bad_line_fails_loud_with_line_number() {
        let input = "chr1\t100\t200\nchr1\tNOT_A_NUMBER\t400\n";
        let err = read_bytes(input.as_bytes()).unwrap_err();
        assert!(err.to_string().contains("line 2"), "{err}");
    }

    #[test]
    fn empty_interval_rejected_at_parse_time() {
        let input = "chr1\t100\t100\n";
        let err = read_bytes(input.as_bytes()).unwrap_err();
        assert!(err.to_string().to_lowercase().contains("empty"), "{err}");
    }

    #[test]
    fn write_bed3_round_trip() {
        let ivs = vec![
            Interval::new("chr1", 100, 200).unwrap(),
            Interval::new("chr1", 300, 400).unwrap(),
        ];
        let mut buf = Vec::new();
        write_bed3(&mut buf, ivs.clone()).unwrap();
        let parsed = read_bytes(&buf).unwrap();
        assert_eq!(parsed, ivs);
    }

    #[test]
    fn write_bed6_includes_strand_and_placeholders() {
        let iv = Interval::with_strand("chr1", 100, 200, Strand::Reverse).unwrap();
        let mut buf = Vec::new();
        write_bed6(&mut buf, [iv]).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert_eq!(s, "chr1\t100\t200\t.\t0\t-\n");
    }

    fn merged(input: &str) -> String {
        let mut out = Vec::new();
        merge_sorted(io::Cursor::new(input.as_bytes()), &mut out).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn merge_sorted_collapses_overlapping_and_touching() {
        let input = "chr1\t100\t200\nchr1\t150\t250\nchr1\t250\t300\nchr1\t400\t500\n\
                     chr2\t10\t20\n";
        assert_eq!(merged(input), "chr1\t100\t300\nchr1\t400\t500\nchr2\t10\t20\n");
    }

    #[test]
    fn merge_sorted_skips_comments_and_blanks() {
        let input = "# h\n\nchr1\t100\t200\nchr1\t150\t250\n";
        assert_eq!(merged(input), "chr1\t100\t250\n");
    }

    #[test]
    fn merge_sorted_rejects_unsorted_same_chrom() {
        let input = "chr1\t100\t200\nchr1\t50\t60\n";
        let mut out = Vec::new();
        let err = merge_sorted(io::Cursor::new(input.as_bytes()), &mut out).unwrap_err();
        assert!(err.to_string().contains("not sorted"), "{err}");
    }

    #[test]
    fn merge_sorted_rejects_inverted_interval() {
        let mut out = Vec::new();
        let err =
            merge_sorted(io::Cursor::new(b"chr1\t200\t100\n".as_slice()), &mut out).unwrap_err();
        assert!(err.to_string().contains("line 1"), "{err}");
    }

    #[test]
    fn merge_sorted_rejects_reappearing_chrom() {
        // chr1 -> chr2 -> chr1: bedtools errors; we must too, not split.
        let input = "chr1\t100\t200\nchr2\t100\t200\nchr1\t300\t400\n";
        let mut out = Vec::new();
        let err = merge_sorted(io::Cursor::new(input.as_bytes()), &mut out).unwrap_err();
        assert!(err.to_string().contains("reappears"), "{err}");
    }

    #[test]
    fn merge_sorted_skips_track_and_browser_preamble() {
        let input = "browser position chr1:1-1000\ntrack name=\"x\"\nchr1\t100\t200\n\
                     chr1\t150\t250\n";
        assert_eq!(merged(input), "chr1\t100\t250\n");
    }

    #[test]
    fn merge_sorted_rejects_coordinate_overflow() {
        let input = "chr1\t0\t184467440737095516160\n"; // > u64::MAX
        let mut out = Vec::new();
        let err = merge_sorted(io::Cursor::new(input.as_bytes()), &mut out).unwrap_err();
        assert!(err.to_string().contains("overflow"), "{err}");
    }

    #[test]
    fn merge_sorted_passes_non_utf8_chrom_through_verbatim() {
        // bedtools treats the contig name as opaque bytes; so do we — no
        // lossy "?" substitution, no spurious error.
        let mut input = b"chr\xff\t100\t200\nchr\xff\t150\t250\n".to_vec();
        let mut out = Vec::new();
        merge_sorted(io::Cursor::new(&mut input), &mut out).unwrap();
        assert_eq!(out, b"chr\xff\t100\t250\n");
    }
}
