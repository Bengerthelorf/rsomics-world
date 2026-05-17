mod detect;
mod parse;
mod reader_gz;

use std::path::Path;

use rsomics_common::Result;

use detect::InputKind;

// id/seq/qual: line terminator (\r?\n) stripped, otherwise verbatim — matches needletail SequenceRecord::id(); parser enforces seq.len() == qual.len()
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OwnedRecord {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

// fail-loud: a parse / IO / truncated-gzip error yields Err then terminates — never silently stops short
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

// detection is by magic bytes only — file extension is ignored
pub fn open_fastq(path: &Path) -> Result<FastqSource> {
    let kind = detect::detect(path)?;
    match kind {
        InputKind::Plain => {
            let r = reader_plain::PlainReader::open(path)?;
            Ok(FastqSource::Plain(r))
        }
        // BGZF = concatenated gzip members with BC subfields; igzip's multi-member path decodes it
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
