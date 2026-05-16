use std::fs::File;
use std::io::Read;
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum InputKind {
    Plain,
    Gz,
    Bgzf,
}

/// Read the first 18 bytes of the file and classify the input format.
///
/// BGZF detection checks for the gzip `BC` extra subfield per the SAM spec
/// §4.1 (BGZF format).  A plain gzip file that happens to start with `0x1f 0x8b`
/// but lacks the BC field is correctly classified as Gz, not Bgzf.  The 18-byte
/// probe window covers the full fixed gzip header (10 bytes) + 2-byte XLEN +
/// the 4-byte BC subfield header (`BC` + 2-byte little-endian size) fully.
pub(crate) fn detect(path: &Path) -> Result<InputKind> {
    let mut f = File::open(path).map_err(|e| {
        RsomicsError::Io(std::io::Error::new(
            e.kind(),
            format!("opening {}: {e}", path.display()),
        ))
    })?;
    // A single `read` may legally return fewer bytes than requested even when
    // more are available; a 1-byte short read here would misclassify a .gz as
    // plain and feed it raw to the parser. Fill the probe across short reads.
    let mut probe = [0u8; 18];
    let mut n = 0;
    while n < probe.len() {
        match f.read(&mut probe[n..]) {
            Ok(0) => break,
            Ok(r) => n += r,
            Err(e) => {
                return Err(RsomicsError::Io(std::io::Error::new(
                    e.kind(),
                    format!("reading header of {}: {e}", path.display()),
                )));
            }
        }
    }

    if n < 2 || probe[0] != 0x1f || probe[1] != 0x8b {
        return Ok(InputKind::Plain);
    }

    // byte 3 = FLG; bit 4 (FEXTRA) must be set.
    // bytes 10–11 = XLEN (little-endian 2 bytes).
    // bytes 12–13 = subfield ID1/ID2 ('B','C').
    // bytes 14–15 = subfield size (must be 2, LE).
    if n >= 16 {
        let flg = probe[3];
        if flg & 0x04 != 0 {
            let xlen = u16::from_le_bytes([probe[10], probe[11]]) as usize;
            if xlen >= 6 && probe[12] == b'B' && probe[13] == b'C' {
                let sf_size = u16::from_le_bytes([probe[14], probe[15]]);
                if sf_size == 2 {
                    return Ok(InputKind::Bgzf);
                }
            }
        }
    }

    Ok(InputKind::Gz)
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    fn write_tmp(bytes: &[u8]) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(bytes).unwrap();
        f.flush().unwrap();
        f
    }

    // Minimal valid BGZF block header:
    // ID1=0x1f ID2=0x8b CM=0x08 FLG=0x04 (FEXTRA)
    // MTIME=0 XFL=0 OS=0xff
    // XLEN=6 (LE) | SI1='B' SI2='C' | SLEN=2 (LE)
    // BSIZE placeholder + padding to make 18 bytes total
    const BGZF_HEADER: &[u8] = &[
        0x1f, 0x8b, // ID1, ID2
        0x08, // CM = deflate
        0x04, // FLG = FEXTRA
        0x00, 0x00, 0x00, 0x00, // MTIME
        0x00, // XFL
        0xff, // OS
        0x06, 0x00, // XLEN = 6
        b'B', b'C', // SI1, SI2
        0x02, 0x00, // SLEN = 2
        0x00, 0x00, // BSIZE placeholder
    ];

    const PLAIN_GZ_HEADER: &[u8] = &[
        0x1f, 0x8b, // ID1, ID2
        0x08, // CM = deflate
        0x00, // FLG = no FEXTRA
        0x00, 0x00, 0x00, 0x00, // MTIME
        0x00, // XFL
        0xff, // OS
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // data padding
    ];

    #[test]
    fn plain_text_detected_as_plain() {
        let f = write_tmp(b"@read1\nACGT\n+\nIIII\n");
        assert_eq!(detect(f.path()).unwrap(), InputKind::Plain);
    }

    #[test]
    fn plain_gz_header_detected_as_gz() {
        let f = write_tmp(PLAIN_GZ_HEADER);
        assert_eq!(detect(f.path()).unwrap(), InputKind::Gz);
    }

    #[test]
    fn bgzf_header_detected_as_bgzf() {
        let f = write_tmp(BGZF_HEADER);
        assert_eq!(detect(f.path()).unwrap(), InputKind::Bgzf);
    }

    #[test]
    fn empty_file_detected_as_plain() {
        let f = write_tmp(b"");
        assert_eq!(detect(f.path()).unwrap(), InputKind::Plain);
    }
}
