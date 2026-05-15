// Gzip permits concatenated members — gunzip/fastp/seqkit/pigz decode them
// transparently. We produce one self-contained gzip member per ~256 KB chunk
// via libdeflate on rayon workers, then write in input order.
// flate2/zlib-rs (~80 MB/s ST) is ~3× slower than fastp's libdeflate path;
// parallel chunking closes that gap.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use libdeflater::{CompressionLvl, Compressor};
use rayon::prelude::*;
use rsomics_common::{Context, Result, RsomicsError};

pub const GZ_CHUNK_BYTES: usize = 256 * 1024;

#[cfg(test)]
pub(crate) const GZ_DEFAULT_LEVEL: i32 = 4;

// 16 × 256 KB = 4 MB peak in-flight. Flush when this fills.
pub const MAX_PENDING_CHUNKS: usize = 16;

/// Compress one buffer to a self-contained gzip member.
fn compress_member(plain: &[u8], level: i32) -> Result<Vec<u8>> {
    let level = CompressionLvl::new(level).map_err(|e| {
        RsomicsError::ConfigError(format!("invalid libdeflate level {level}: {e:?}"))
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
pub fn write_chunks_gz<W: Write>(out: &mut W, chunks: Vec<Vec<u8>>, level: i32) -> Result<()> {
    let compressed: Vec<Result<Vec<u8>>> = chunks
        .into_par_iter()
        .map(|c| compress_member(&c, level))
        .collect();
    for c in compressed {
        let bytes = c?;
        out.write_all(&bytes).rs_context("writing gzip member")?;
    }
    Ok(())
}

/// Write one record in `@id\nseq\n+\nqual\n` form. `id` carries no leading `@`.
fn write_plain_fastq_record<W: Write>(
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

/// Append-style writer that buffers plain bytes until a chunk fills, then
/// emits the chunk via the parallel-gz pipeline. Used by `pipeline.rs`'s
/// SE/PE write paths. Wraps a `BufWriter` so plain-text output stays
/// fast too.
pub struct ChunkedWriter {
    inner: BufWriter<File>,
    buffer: Vec<u8>,
    gzipped: bool,
    pending_chunks: Vec<Vec<u8>>,
    level: i32,
}

impl ChunkedWriter {
    /// Open `path` for writing. `.gz` extension selects parallel-gz at
    /// `level`; any other extension writes plain bytes and ignores
    /// `level`.
    ///
    /// # Errors
    ///
    /// `Io` if the file cannot be created.
    pub fn create(path: &Path, level: i32) -> Result<Self> {
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
            level,
        })
    }

    /// Splits off a chunk and queues it for compression when the buffer
    /// crosses [`GZ_CHUNK_BYTES`]; if more than [`MAX_PENDING_CHUNKS`]
    /// chunks are waiting, flushes the queue so the writer's footprint
    /// stays bounded regardless of input size.
    ///
    /// # Errors
    ///
    /// `Io` if a plain-text write fails; `UpstreamError` if libdeflate
    /// compression fails while the pending queue is being drained.
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
                if self.pending_chunks.len() >= MAX_PENDING_CHUNKS {
                    self.drain_pending()?;
                }
            }
        } else {
            write_plain_fastq_record(&mut self.inner, id, seq, qual)
                .rs_context("writing plain FASTQ record")?;
        }
        Ok(())
    }

    /// Compress the currently queued chunks in parallel and write them
    /// to disk, leaving the in-progress `buffer` untouched. Called both
    /// during a run (to bound memory) and from `finalize` for the tail.
    fn drain_pending(&mut self) -> Result<()> {
        if self.pending_chunks.is_empty() {
            return Ok(());
        }
        let chunks = std::mem::take(&mut self.pending_chunks);
        write_chunks_gz(&mut self.inner, chunks, self.level)
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
        self.drain_pending()?;
        self.inner.flush().rs_context("flushing gz writer")?;
        Ok(())
    }

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
        let mut w = ChunkedWriter::create(tmp.path(), GZ_DEFAULT_LEVEL).unwrap();
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
        let mut w = ChunkedWriter::create(tmp.path(), GZ_DEFAULT_LEVEL).unwrap();
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
