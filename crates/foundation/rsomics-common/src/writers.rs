//! Output sink helpers that abstract "write to a file path OR to stdout".
//!
//! BAM/CRAM/SAM tools idiomatically pipe — `samtools view in.bam | samtools
//! sort -o out.bam -` — so every subcommand that produces a binary stream
//! needs the same `Option<&Path>` dispatch. Centralising it here keeps each
//! tool's I/O boundary one line at the sink site.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::error::{Context, Result};

/// Open an output writer: `Some(path)` → buffered file; `None` → buffered
/// stdout. The returned writer is `Box<dyn Write + Send>` so callers can
/// pass it through trait-object boundaries (parallel scatter-gather
/// pipelines, generic record writers).
///
/// # Errors
///
/// Returns `RsomicsError::Io` if `path` is `Some(p)` and `p` cannot be
/// created or truncated.
pub fn open_output(path: Option<&Path>) -> Result<Box<dyn Write + Send>> {
    match path {
        Some(p) => {
            let f =
                File::create(p).rs_with_context(|| format!("creating output {}", p.display()))?;
            Ok(Box::new(BufWriter::new(f)))
        }
        None => Ok(Box::new(BufWriter::new(io::stdout()))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    #[test]
    fn writes_to_file_when_path_supplied() {
        let tmp = tempfile::NamedTempFile::new().expect("tempfile");
        {
            let mut w = open_output(Some(tmp.path())).expect("open");
            w.write_all(b"hello\n").expect("write");
            w.flush().expect("flush");
        }
        let mut buf = String::new();
        std::fs::File::open(tmp.path())
            .expect("open back")
            .read_to_string(&mut buf)
            .expect("read");
        assert_eq!(buf, "hello\n");
    }

    #[test]
    fn returns_stdout_writer_when_path_is_none() {
        // Verify the call succeeds and returns a usable writer — we don't
        // capture stdout in unit tests, but the buffer-write itself must
        // not error.
        let mut w = open_output(None).expect("open stdout");
        w.write_all(b"").expect("empty write OK");
    }
}
