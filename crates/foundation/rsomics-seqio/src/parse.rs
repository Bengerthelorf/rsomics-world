use std::io::BufRead;

use rsomics_common::{Result, RsomicsError};

use crate::OwnedRecord;

/// Parse one FASTQ record from `reader`.
///
/// Returns `Ok(None)` at clean EOF (no bytes consumed for this record).
/// Returns `Err` for truncated input, missing `@`/`+` markers, or seq/qual
/// length mismatch.  Trailing `\r` on any line is stripped; otherwise line
/// content is kept verbatim.
pub(crate) fn parse_record<R: BufRead>(reader: &mut R) -> Result<Option<OwnedRecord>> {
    let mut header = Vec::new();
    let n = reader
        .read_until(b'\n', &mut header)
        .map_err(RsomicsError::Io)?;
    if n == 0 {
        return Ok(None);
    }
    trim_newline(&mut header);
    if header.first() != Some(&b'@') {
        let got = match header.first() {
            None => "empty line".to_string(),
            Some(&b) => format!("{:?}", b as char),
        };
        return Err(RsomicsError::InvalidInput(format!(
            "expected '@' at record start, got {got}"
        )));
    }
    let id = header[1..].to_vec();

    let mut seq = Vec::new();
    if reader
        .read_until(b'\n', &mut seq)
        .map_err(RsomicsError::Io)?
        == 0
    {
        return Err(RsomicsError::InvalidInput(
            "truncated FASTQ: missing sequence line".into(),
        ));
    }
    trim_newline(&mut seq);

    let mut sep = Vec::new();
    if reader
        .read_until(b'\n', &mut sep)
        .map_err(RsomicsError::Io)?
        == 0
    {
        return Err(RsomicsError::InvalidInput(
            "truncated FASTQ: missing '+' separator".into(),
        ));
    }
    trim_newline(&mut sep);
    if sep.first() != Some(&b'+') {
        return Err(RsomicsError::InvalidInput(format!(
            "expected '+' separator, got {:?}",
            sep.first().copied().unwrap_or(0) as char
        )));
    }

    let mut qual = Vec::new();
    if reader
        .read_until(b'\n', &mut qual)
        .map_err(RsomicsError::Io)?
        == 0
    {
        return Err(RsomicsError::InvalidInput(
            "truncated FASTQ: missing quality line".into(),
        ));
    }
    trim_newline(&mut qual);

    if seq.len() != qual.len() {
        return Err(RsomicsError::InvalidInput(format!(
            "seq/qual length mismatch: {} vs {}",
            seq.len(),
            qual.len()
        )));
    }

    Ok(Some(OwnedRecord { id, seq, qual }))
}

#[inline]
fn trim_newline(buf: &mut Vec<u8>) {
    if buf.last() == Some(&b'\n') {
        buf.pop();
    }
    if buf.last() == Some(&b'\r') {
        buf.pop();
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    fn parse(bytes: &[u8]) -> Result<Option<OwnedRecord>> {
        parse_record(&mut Cursor::new(bytes))
    }

    #[test]
    fn basic_record() {
        let rec = parse(b"@read1\nACGT\n+\nIIII\n").unwrap().unwrap();
        assert_eq!(rec.id, b"read1");
        assert_eq!(rec.seq, b"ACGT");
        assert_eq!(rec.qual, b"IIII");
    }

    #[test]
    fn crlf_stripped() {
        let rec = parse(b"@read1\r\nACGT\r\n+\r\nIIII\r\n").unwrap().unwrap();
        assert_eq!(rec.id, b"read1");
        assert_eq!(rec.seq, b"ACGT");
        assert_eq!(rec.qual, b"IIII");
    }

    #[test]
    fn header_with_spaces_preserved() {
        let rec = parse(b"@read1 desc here\nACGT\n+\nIIII\n")
            .unwrap()
            .unwrap();
        assert_eq!(rec.id, b"read1 desc here");
    }

    #[test]
    fn empty_input_returns_none() {
        assert!(parse(b"").unwrap().is_none());
    }

    #[test]
    fn empty_record() {
        let rec = parse(b"@\n\n+\n\n").unwrap().unwrap();
        assert_eq!(rec.id, b"");
        assert_eq!(rec.seq, b"");
        assert_eq!(rec.qual, b"");
    }

    #[test]
    fn seq_qual_mismatch_is_error() {
        let e = parse(b"@r\nACGT\n+\nIII\n").unwrap_err();
        assert!(matches!(e, RsomicsError::InvalidInput(_)));
    }

    #[test]
    fn missing_at_is_error() {
        let e = parse(b"read1\nACGT\n+\nIIII\n").unwrap_err();
        assert!(matches!(e, RsomicsError::InvalidInput(_)));
    }

    #[test]
    fn truncated_after_header_is_error() {
        let e = parse(b"@read1\n").unwrap_err();
        assert!(matches!(e, RsomicsError::InvalidInput(_)));
    }

    #[test]
    fn missing_plus_separator_is_error() {
        let e = parse(b"@r\nACGT\n-\nIIII\n").unwrap_err();
        assert!(matches!(e, RsomicsError::InvalidInput(_)));
    }
}
