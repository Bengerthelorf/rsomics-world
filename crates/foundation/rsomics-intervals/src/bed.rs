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
}
