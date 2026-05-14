//! Parallel-gzip output via rayon-distributed libdeflate compression.
//!
//! `flate2` with the `zlib-rs` backend gives us ~80 MB/s single-threaded
//! compression — competitive with system gzip but ~3× slower than fastp's
//! libdeflate-backed write path. Compression dominates the wall clock for
//! a no-trim pass (decompress takes ~13 ms on a 17 MB input; compress
//! 380 ms), so this is the only axis where a parallel-codec win is large
//! enough to matter.
//!
//! ## Approach
//!
//! Gzip permits concatenated members in one file — `gunzip` / fastp /
//! seqkit / pigz all decode it transparently. So we batch trimmed FASTQ
//! bytes into ~256 KB chunks, hand each chunk to a rayon worker that
//! produces a self-contained gzip member via libdeflate, then write the
//! compressed members to the output file in input order.
//!
//! For plain-text output we skip this whole module and write directly
//! through a `BufWriter`.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use libdeflater::{CompressionLvl, Compressor};
use rayon::prelude::*;
use rsomics_common::{Context, Result, RsomicsError};

/// Target plain-FASTQ bytes per compression chunk. 256 KB gives ~4×
/// gzip ratio → ~64 KB compressed members, modest header overhead, and
/// enough work per rayon dispatch to amortise the per-call cost.
pub const GZ_CHUNK_BYTES: usize = 256 * 1024;

/// libdeflate compression level. fastp defaults to level 4 (out of 12)
/// for a good ratio/speed trade-off. We match.
pub const GZ_LEVEL: i32 = 4;

/// Compress one buffer to a self-contained gzip member.
fn compress_member(plain: &[u8]) -> Result<Vec<u8>> {
    let level = CompressionLvl::new(GZ_LEVEL).map_err(|e| {
        RsomicsError::ConfigError(format!("invalid libdeflate level {GZ_LEVEL}: {e:?}"))
    })?;
    let mut compressor = Compressor::new(level);
    let bound = compressor.gzip_compress_bound(plain.len());
    let mut out = vec![0u8; bound];
    let n = compressor
        .gzip_compress(plain, &mut out)
        .map_err(|e| RsomicsError::UpstreamError(format!("libdeflate gzip_compress: {e:?}")))?;
    out.truncate(n);
    Ok(out)
}

/// Compress a list of plain-byte chunks in parallel and write the
/// resulting gzip members to `out` in input order.
///
/// # Errors
///
/// `UpstreamError` if libdeflate compression fails; `Io` if the write to
/// `out` fails.
pub fn write_chunks_gz<W: Write>(out: &mut W, chunks: Vec<Vec<u8>>) -> Result<()> {
    let compressed: Vec<Result<Vec<u8>>> = chunks
        .into_par_iter()
        .map(|c| compress_member(&c))
        .collect();
    for c in compressed {
        let bytes = c?;
        out.write_all(&bytes).rs_context("writing gzip member")?;
    }
    Ok(())
}

/// Append-style writer that buffers plain bytes until a chunk fills, then
/// emits the chunk via the parallel-gz pipeline. Used by `pipeline.rs`'s
/// SE/PE write paths. Wraps a `BufWriter` so plain-text output stays
/// fast too.
pub struct ChunkedWriter {
    inner: BufWriter<File>,
    buffer: Vec<u8>,
    gzipped: bool,
    pending_chunks: Vec<Vec<u8>>,
}

impl ChunkedWriter {
    /// Open `path` for writing. `.gz` extension selects parallel-gz; any
    /// other extension writes plain bytes through a `BufWriter`.
    ///
    /// # Errors
    ///
    /// `Io` if the file cannot be created.
    pub fn create(path: &Path) -> Result<Self> {
        let f =
            File::create(path).rs_with_context(|| format!("creating output {}", path.display()))?;
        let gzipped = path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("gz"));
        Ok(Self {
            inner: BufWriter::with_capacity(GZ_CHUNK_BYTES * 2, f),
            buffer: Vec::with_capacity(GZ_CHUNK_BYTES + 16 * 1024),
            gzipped,
            pending_chunks: Vec::new(),
        })
    }

    /// Append a complete plain-byte record `@id\nseq\n+\nqual\n` to the
    /// buffer. Splits off a chunk and queues it for compression when the
    /// buffer crosses [`GZ_CHUNK_BYTES`].
    ///
    /// # Errors
    ///
    /// `Io` if a plain-text write fails.
    pub fn write_record(&mut self, id: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
        if self.gzipped {
            self.buffer.push(b'@');
            self.buffer.extend_from_slice(id);
            self.buffer.push(b'\n');
            self.buffer.extend_from_slice(seq);
            self.buffer.extend_from_slice(b"\n+\n");
            self.buffer.extend_from_slice(qual);
            self.buffer.push(b'\n');
            if self.buffer.len() >= GZ_CHUNK_BYTES {
                let full = std::mem::take(&mut self.buffer);
                self.buffer.reserve(GZ_CHUNK_BYTES + 16 * 1024);
                self.pending_chunks.push(full);
            }
        } else {
            self.inner.write_all(b"@").rs_context("write")?;
            self.inner.write_all(id).rs_context("write")?;
            self.inner.write_all(b"\n").rs_context("write")?;
            self.inner.write_all(seq).rs_context("write")?;
            self.inner.write_all(b"\n+\n").rs_context("write")?;
            self.inner.write_all(qual).rs_context("write")?;
            self.inner.write_all(b"\n").rs_context("write")?;
        }
        Ok(())
    }

    /// Compress all pending chunks and flush to disk. Idempotent on
    /// subsequent calls.
    ///
    /// # Errors
    ///
    /// `Io` or `UpstreamError` if compression / write fails.
    pub fn flush_pending(&mut self) -> Result<()> {
        if !self.gzipped {
            self.inner.flush().rs_context("flushing plain writer")?;
            return Ok(());
        }
        if !self.buffer.is_empty() {
            let full = std::mem::take(&mut self.buffer);
            self.pending_chunks.push(full);
        }
        let chunks = std::mem::take(&mut self.pending_chunks);
        write_chunks_gz(&mut self.inner, chunks)?;
        self.inner.flush().rs_context("flushing gz writer")?;
        Ok(())
    }

    /// Final flush + close. Drops the writer.
    ///
    /// # Errors
    ///
    /// Same as [`Self::flush_pending`].
    pub fn finalize(mut self) -> Result<()> {
        self.flush_pending()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    fn read_all(p: &Path) -> Vec<u8> {
        let mut bytes = Vec::new();
        File::open(p).unwrap().read_to_end(&mut bytes).unwrap();
        bytes
    }

    #[test]
    fn plain_round_trips() {
        let tmp = tempfile::Builder::new().suffix(".fq").tempfile().unwrap();
        let mut w = ChunkedWriter::create(tmp.path()).unwrap();
        w.write_record(b"r1", b"ACGT", b"IIII").unwrap();
        w.finalize().unwrap();
        assert_eq!(read_all(tmp.path()), b"@r1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn gzipped_output_starts_with_gzip_magic() {
        let tmp = tempfile::Builder::new()
            .suffix(".fq.gz")
            .tempfile()
            .unwrap();
        let mut w = ChunkedWriter::create(tmp.path()).unwrap();
        for _ in 0..1000 {
            w.write_record(b"r", b"ACGTACGTACGTACGTACGT", b"IIIIIIIIIIIIIIIIIIII")
                .unwrap();
        }
        w.finalize().unwrap();
        let bytes = read_all(tmp.path());
        assert_eq!(&bytes[..2], &[0x1f, 0x8b], "gzip magic bytes");
        let mut gz = flate2::read::MultiGzDecoder::new(&bytes[..]);
        let mut plain = Vec::new();
        gz.read_to_end(&mut plain).unwrap();
        assert_eq!(plain.len(), 1000 * (1 + 1 + 1 + 20 + 3 + 20 + 1));
    }
}
