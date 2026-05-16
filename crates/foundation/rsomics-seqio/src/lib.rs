mod detect;
mod parse;
mod reader_gz;

use std::path::Path;

use rsomics_common::Result;

use detect::InputKind;

/// One FASTQ record, fully owned.
///
/// `id` is the header line after `@` with the trailing `\n` and a preceding
/// `\r` (if any) stripped, otherwise verbatim — spaces, slashes, and
/// description fields are preserved (matching needletail's
/// `SequenceRecord::id()`). `seq` and `qual` are the sequence and quality
/// line bytes with the same line-terminator stripping, otherwise verbatim.
/// The parser enforces `seq.len() == qual.len()`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OwnedRecord {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

/// An iterator of [`OwnedRecord`]s over a FASTQ source.
///
/// Created by [`open_fastq`].  On any parse or I/O error the iterator yields
/// `Err(RsomicsError)` and then terminates.  Truncated or corrupt gzip data
/// surfaces as `Err` rather than silently stopping.
pub enum FastqSource {
    Plain(reader_plain::PlainReader),
    Gz(reader_gz::GzReader),
}

impl Iterator for FastqSource {
    type Item = Result<OwnedRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            FastqSource::Plain(r) => r.next(),
            FastqSource::Gz(r) => r.next(),
        }
    }
}

/// Open a FASTQ file and return an iterator of owned records.
///
/// Detection is by magic bytes only — the file extension is ignored:
/// - gzip (`0x1f 0x8b`), including BGZF (its `BC`-subfield blocks are
///   concatenated gzip members the igzip backend decodes transparently):
///   `FastqSource::Gz`, a dedicated decode thread + parallel slab parse.
/// - Anything else: `FastqSource::Plain` backed by a `BufReader`.
///
/// # Errors
///
/// Returns `RsomicsError::Io` if the file cannot be opened or if the first
/// bytes cannot be read for detection.
pub fn open_fastq(path: &Path) -> Result<FastqSource> {
    let kind = detect::detect(path)?;
    match kind {
        InputKind::Plain => {
            let r = reader_plain::PlainReader::open(path)?;
            Ok(FastqSource::Plain(r))
        }
        // BGZF is a stream of concatenated gzip members with BC extra subfields;
        // the igzip backend's multi-member handling decodes it correctly.
        InputKind::Gz | InputKind::Bgzf => {
            let r = reader_gz::GzReader::open(path)?;
            Ok(FastqSource::Gz(r))
        }
    }
}

pub(crate) mod reader_plain {
    use std::fs::File;
    use std::io::BufReader;
    use std::path::Path;

    use rsomics_common::{Result, RsomicsError};

    use crate::OwnedRecord;
    use crate::parse::parse_record;

    const BUF: usize = 256 * 1024;

    pub struct PlainReader {
        inner: BufReader<File>,
        done: bool,
    }

    impl PlainReader {
        pub fn open(path: &Path) -> Result<Self> {
            let f = File::open(path).map_err(|e| {
                RsomicsError::Io(std::io::Error::new(
                    e.kind(),
                    format!("opening {}: {e}", path.display()),
                ))
            })?;
            Ok(Self {
                inner: BufReader::with_capacity(BUF, f),
                done: false,
            })
        }
    }

    impl Iterator for PlainReader {
        type Item = Result<OwnedRecord>;

        fn next(&mut self) -> Option<Self::Item> {
            if self.done {
                return None;
            }
            match parse_record(&mut self.inner) {
                Ok(None) => {
                    self.done = true;
                    None
                }
                Ok(Some(r)) => Some(Ok(r)),
                Err(e) => {
                    self.done = true;
                    Some(Err(e))
                }
            }
        }
    }
}
