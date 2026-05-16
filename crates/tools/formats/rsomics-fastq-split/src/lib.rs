use std::path::{Path, PathBuf};

use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::{OwnedRecord, open_fastq};
use serde::Serialize;

use rsomics_fqgz::ChunkedWriter;

/// How to partition the input.
///
/// fastp 0.20.1 exposes `--split N` (byFileNumber) and `--split_by_lines L`
/// (byFileLines). fastp's byFileNumber distributes by a *file-size estimate*
/// (not an exact read count), so it is not a deterministic byte-compat oracle;
/// `ByNumber` here is **exact-count** (one counting pass → `ceil(total/N)`
/// contiguous reads per file), a deliberate, documented improvement. `ByLines`
/// is deterministic on both sides and is byte-equal to fastp 0.20.1.
#[derive(Debug, Clone, Copy)]
pub enum SplitMode {
    /// `--split_by_lines L`: each file gets `L/4` reads (last file the
    /// remainder). `L` must be a multiple of 4 (a FASTQ record is 4 lines).
    ByLines(usize),
    /// `--split N`: exactly enough files of `ceil(total/N)` reads each.
    ByNumber(usize),
}

#[derive(Debug, Clone)]
pub struct SplitConfig {
    pub mode: SplitMode,
    /// fastp `--split_prefix_digits` (zero-pad width of the numeric prefix).
    pub digits: usize,
}

/// fastp 0.20.1 split-file name: `<zero-padded 1-based index>.<basename(out)>`
/// in `out`'s parent directory (so `out.fq.gz` → `0001.out.fq.gz`, keeping the
/// `.gz` suffix so the per-file writer still gzip-compresses).
fn split_path(out: &Path, idx: usize, digits: usize) -> PathBuf {
    let base = out.file_name().map_or_else(
        || std::ffi::OsString::from("out"),
        std::ffi::OsStr::to_os_string,
    );
    let name = format!("{idx:0width$}.{}", base.to_string_lossy(), width = digits);
    match out.parent() {
        Some(p) if !p.as_os_str().is_empty() => p.join(name),
        _ => PathBuf::from(name),
    }
}

struct RollingWriter<'a> {
    out: &'a Path,
    digits: usize,
    compression: i32,
    per_file: usize,
    idx: usize,
    in_file: usize,
    cur: Option<ChunkedWriter>,
    files_written: usize,
}

impl<'a> RollingWriter<'a> {
    fn new(out: &'a Path, digits: usize, compression: i32, per_file: usize) -> Self {
        Self {
            out,
            digits,
            compression,
            per_file,
            idx: 1,
            in_file: 0,
            cur: None,
            files_written: 0,
        }
    }

    fn ensure_open(&mut self) -> Result<()> {
        if self.cur.is_none() {
            let path = split_path(self.out, self.idx, self.digits);
            self.cur = Some(ChunkedWriter::create(&path, self.compression)?);
            self.files_written += 1;
        }
        Ok(())
    }

    fn write(&mut self, rec: &OwnedRecord) -> Result<()> {
        if self.in_file == self.per_file {
            self.roll()?;
        }
        self.ensure_open()?;
        self.cur
            .as_mut()
            .expect("writer opened by ensure_open")
            .write_record(&rec.id, &rec.seq, &rec.qual)?;
        self.in_file += 1;
        Ok(())
    }

    fn roll(&mut self) -> Result<()> {
        if let Some(w) = self.cur.take() {
            w.finalize()?;
        }
        self.idx += 1;
        self.in_file = 0;
        Ok(())
    }

    fn finalize(mut self) -> Result<usize> {
        if let Some(w) = self.cur.take() {
            w.finalize()?;
        }
        Ok(self.files_written)
    }
}

#[derive(Debug, Default, Clone, Serialize)]
pub struct SplitReport {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mode: Option<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub input_r1: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub input_r2: Option<String>,
    pub reads_in: u64,
    pub bases_in: u64,
    pub files_written: u64,
}

fn reads_per_file(mode: SplitMode, total: Option<usize>) -> Result<usize> {
    match mode {
        SplitMode::ByLines(l) => {
            if l == 0 || l % 4 != 0 {
                return Err(RsomicsError::ConfigError(format!(
                    "--split_by_lines must be a positive multiple of 4 (a FASTQ record is 4 lines), got {l}"
                )));
            }
            Ok(l / 4)
        }
        SplitMode::ByNumber(n) => {
            if n == 0 {
                return Err(RsomicsError::ConfigError("--split must be > 0".into()));
            }
            let total = total.expect("ByNumber requires a counted total");
            Ok(total.div_ceil(n).max(1))
        }
    }
}

fn count_reads(input: &Path) -> Result<usize> {
    let mut n = 0usize;
    for r in open_fastq(input)? {
        r?;
        n += 1;
    }
    Ok(n)
}

pub struct Pipeline<'cfg> {
    pub cfg: &'cfg SplitConfig,
    pub compression: i32,
}

impl<'cfg> Pipeline<'cfg> {
    #[must_use]
    pub fn new(cfg: &'cfg SplitConfig, compression: i32) -> Self {
        Self { cfg, compression }
    }

    /// # Errors
    ///
    /// Propagates input parse / output write errors and config errors.
    pub fn run_se(&self, input: &Path, out: &Path) -> Result<SplitReport> {
        let total = match self.cfg.mode {
            SplitMode::ByNumber(_) => Some(count_reads(input)?),
            SplitMode::ByLines(_) => None,
        };
        let per_file = reads_per_file(self.cfg.mode, total)?;
        let mut w = RollingWriter::new(out, self.cfg.digits, self.compression, per_file);
        let mut report = SplitReport {
            mode: Some("SE"),
            input_r1: Some(input.display().to_string()),
            ..SplitReport::default()
        };
        for r in open_fastq(input)? {
            let rec = r?;
            report.reads_in += 1;
            report.bases_in += rec.seq.len() as u64;
            w.write(&rec)?;
        }
        report.files_written = w.finalize()? as u64;
        Ok(report)
    }

    /// # Errors
    ///
    /// Propagates parse / write / config errors; errors if the two inputs have
    /// a differing record count.
    pub fn run_pe(&self, in1: &Path, in2: &Path, out1: &Path, out2: &Path) -> Result<SplitReport> {
        let total = match self.cfg.mode {
            SplitMode::ByNumber(_) => Some(count_reads(in1)?),
            SplitMode::ByLines(_) => None,
        };
        let per_file = reads_per_file(self.cfg.mode, total)?;
        let mut w1 = RollingWriter::new(out1, self.cfg.digits, self.compression, per_file);
        let mut w2 = RollingWriter::new(out2, self.cfg.digits, self.compression, per_file);
        let mut report = SplitReport {
            mode: Some("PE"),
            input_r1: Some(in1.display().to_string()),
            input_r2: Some(in2.display().to_string()),
            ..SplitReport::default()
        };
        let mut r1 = open_fastq(in1)?;
        let mut r2 = open_fastq(in2)?;
        loop {
            match (r1.next(), r2.next()) {
                (Some(a), Some(b)) => {
                    let (ra, rb) = (a?, b?);
                    report.reads_in += 2;
                    report.bases_in += (ra.seq.len() + rb.seq.len()) as u64;
                    w1.write(&ra)?;
                    w2.write(&rb)?;
                }
                (None, None) => break,
                _ => {
                    return Err(RsomicsError::InvalidInput(
                        "PE input record counts diverge".into(),
                    ));
                }
            }
        }
        let f1 = w1.finalize()?;
        let f2 = w2.finalize()?;
        debug_assert_eq!(f1, f2, "PE split must produce matching file counts");
        report.files_written = (f1 + f2) as u64;
        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_path_zero_pads_and_keeps_dir_and_gz_suffix() {
        let p = split_path(Path::new("/tmp/x/out.fq.gz"), 1, 4);
        assert_eq!(p, Path::new("/tmp/x/0001.out.fq.gz"));
        let p2 = split_path(Path::new("out.fq"), 12, 4);
        assert_eq!(p2, Path::new("0012.out.fq"));
        let p3 = split_path(Path::new("/tmp/o.fq"), 7, 3);
        assert_eq!(p3, Path::new("/tmp/007.o.fq"));
    }

    #[test]
    fn by_lines_must_be_multiple_of_four() {
        assert!(reads_per_file(SplitMode::ByLines(10), None).is_err());
        assert!(reads_per_file(SplitMode::ByLines(0), None).is_err());
        assert_eq!(reads_per_file(SplitMode::ByLines(8), None).unwrap(), 2);
        assert_eq!(reads_per_file(SplitMode::ByLines(400), None).unwrap(), 100);
    }

    #[test]
    fn by_number_is_exact_ceil() {
        assert_eq!(reads_per_file(SplitMode::ByNumber(4), Some(10)).unwrap(), 3);
        assert_eq!(reads_per_file(SplitMode::ByNumber(3), Some(9)).unwrap(), 3);
        assert_eq!(reads_per_file(SplitMode::ByNumber(10), Some(3)).unwrap(), 1);
        assert_eq!(reads_per_file(SplitMode::ByNumber(5), Some(0)).unwrap(), 1);
        assert!(reads_per_file(SplitMode::ByNumber(0), Some(10)).is_err());
    }
}
