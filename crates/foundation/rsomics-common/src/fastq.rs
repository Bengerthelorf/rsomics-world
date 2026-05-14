//! Shared FASTQ output writer and per-record formatter.
//!
//! Every FASTQ-producing crate in the workspace (`rsomics-fastq-trim`,
//! `rsomics-fastq-quality`, `rsomics-fastq-umi`, …) writes the same
//! `@id\nseq\n+\nqual\n` format and needs a gzip-aware sink. Putting both
//! in one module avoids three near-identical duplicates and keeps the
//! gzip-finalisation contract in one place.
//!
//! ## Why an enum, not a `Box<dyn Write>`
//!
//! The plain path is the common one. Enum dispatch keeps it monomorphic
//! and lets the buffered writer's small writes inline through `BufWriter`
//! without virtual-call overhead per record.
//!
//! ## Finalisation contract
//!
//! `GzEncoder::Drop` calls `try_finish` and swallows any late error (disk
//! full at the trailer write). Callers should prefer the explicit
//! [`Writer::finalize`] consuming-method so trailer-write errors surface
//! as `RsomicsError::Io`.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use flate2::Compression;
use flate2::write::GzEncoder;

use crate::error::{Context, Result};

/// Gzip-aware FASTQ output writer. Auto-selects gzipped output when the
/// path ends in `.gz`; plain otherwise. Both arms wrap a `BufWriter` so
/// needletail's small per-record writes batch into larger I/O.
pub enum Writer {
    Plain(BufWriter<File>),
    // Boxed so the enum's variant-size disparity (a GzEncoder is much
    // bigger than a BufWriter alone) doesn't bloat every `Plain` value
    // by carrying the gzip-state space.
    Gzip(Box<GzEncoder<BufWriter<File>>>),
}

impl Writer {
    /// Open `path` for writing. Gzip is selected by `.gz` extension match
    /// (case-insensitive).
    ///
    /// # Errors
    ///
    /// Returns `RsomicsError::Io` if the file cannot be created.
    pub fn create(path: &Path) -> Result<Self> {
        let file = File::create(path)
            .rs_with_context(|| format!("creating output FASTQ {}", path.display()))?;
        let buf = BufWriter::new(file);
        if path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("gz"))
        {
            Ok(Self::Gzip(Box::new(GzEncoder::new(buf, Compression::default()))))
        } else {
            Ok(Self::Plain(buf))
        }
    }

    /// Flush + finish. For the gzip variant this writes the trailer; the
    /// `Drop` impl would silently swallow trailer-write errors so callers
    /// that care about durable output should always end with this.
    ///
    /// # Errors
    ///
    /// Returns `RsomicsError::Io` if the underlying flush / finish fails.
    pub fn finalize(self) -> Result<()> {
        match self {
            Self::Plain(mut w) => w.flush().rs_context("flushing plain output writer")?,
            Self::Gzip(w) => {
                w.finish().rs_context("finishing gzip output stream")?;
            }
        }
        Ok(())
    }
}

impl Write for Writer {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Plain(w) => w.write(buf),
            Self::Gzip(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Plain(w) => w.flush(),
            Self::Gzip(w) => w.flush(),
        }
    }
}

/// Write one FASTQ record to `w` in canonical `@id\nseq\n+\nqual\n` form.
/// `id` must NOT include the leading `@`; this fn adds it.
///
/// # Errors
///
/// Forwards any `std::io::Error` from the underlying writer.
pub fn write_record<W: Write>(
    w: &mut W,
    id: &[u8],
    seq: &[u8],
    qual: &[u8],
) -> std::io::Result<()> {
    w.write_all(b"@")?;
    w.write_all(id)?;
    w.write_all(b"\n")?;
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    w.write_all(qual)?;
    w.write_all(b"\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    #[test]
    fn plain_writer_round_trips() {
        let tmp = tempfile::Builder::new().suffix(".fq").tempfile().unwrap();
        let path = tmp.path().to_path_buf();
        let mut w = Writer::create(&path).unwrap();
        write_record(&mut w, b"r1", b"ACGT", b"IIII").unwrap();
        w.finalize().unwrap();
        let mut content = String::new();
        File::open(&path)
            .unwrap()
            .read_to_string(&mut content)
            .unwrap();
        assert_eq!(content, "@r1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn gzip_writer_writes_magic_bytes() {
        let tmp = tempfile::Builder::new()
            .suffix(".fq.gz")
            .tempfile()
            .unwrap();
        let path = tmp.path().to_path_buf();
        let mut w = Writer::create(&path).unwrap();
        write_record(&mut w, b"r1", b"ACGT", b"IIII").unwrap();
        w.finalize().unwrap();
        let mut bytes = Vec::new();
        File::open(&path).unwrap().read_to_end(&mut bytes).unwrap();
        assert_eq!(&bytes[..2], &[0x1f, 0x8b], "gzip magic bytes");
    }
}
