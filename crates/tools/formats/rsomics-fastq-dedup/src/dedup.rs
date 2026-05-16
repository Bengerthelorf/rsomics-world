use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use flate2::Compression;
use flate2::write::GzEncoder;
use needletail::parse_fastx_file;
use rsomics_common::{Context, Result, RsomicsError};
use rsomics_kmer::{KmerIter, murmur3_x64_128};

use crate::report::DedupReport;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DedupMode {
    KmerBin,
    FullSeq,
}

#[derive(Debug, Clone, Copy)]
pub struct DedupConfig {
    pub mode: DedupMode,
    pub k: usize,
    pub tail_offset: usize,
}

impl Default for DedupConfig {
    fn default() -> Self {
        Self {
            mode: DedupMode::KmerBin,
            k: 12,
            tail_offset: 0,
        }
    }
}

pub fn run_se(input: &Path, output: &Path, cfg: DedupConfig) -> Result<DedupReport> {
    let mut reader = parse_fastx_file(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("opening {}: {e}", input.display())))?;
    let mut writer = open_writer(output)?;

    let mut seen: HashSet<u128> = HashSet::new();
    let mut report = DedupReport {
        input: Some(input.display().to_string()),
        output: Some(output.display().to_string()),
        ..DedupReport::default()
    };

    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| RsomicsError::InvalidInput(format!("FASTQ parse: {e}")))?;
        report.reads_in += 1;
        let key = match cfg.mode {
            DedupMode::FullSeq => full_seq_key(&rec.seq()),
            DedupMode::KmerBin => kmer_bin_key(&rec.seq(), cfg.k, cfg.tail_offset)?,
        };
        if seen.insert(key) {
            let qual = rec
                .qual()
                .ok_or_else(|| RsomicsError::InvalidInput("FASTQ record missing quality".into()))?;
            write_record(&mut writer, rec.id(), &rec.seq(), qual)?;
            report.reads_out += 1;
        } else {
            report.duplicates_removed += 1;
        }
    }
    writer.flush().rs_context("flushing output")?;
    Ok(report)
}

fn full_seq_key(seq: &[u8]) -> u128 {
    murmur3_x64_128(seq, 0)
}

fn kmer_bin_key(seq: &[u8], k: usize, tail_offset: usize) -> Result<u128> {
    if seq.len() < k {
        return Ok(murmur3_x64_128(seq, 0));
    }
    let start = seq.len().saturating_sub(k + tail_offset);
    let window = &seq[start..start + k];
    let mut it = KmerIter::new(window, k, true)
        .map_err(|e| RsomicsError::InvalidInput(format!("k-mer window: {e}")))?;
    match it.next() {
        Some(Ok(kmer)) => Ok(u128::from(kmer)),
        _ => Ok(murmur3_x64_128(seq, 0)),
    }
}

fn open_writer(path: &Path) -> Result<Box<dyn Write>> {
    let file =
        File::create(path).rs_with_context(|| format!("creating output {}", path.display()))?;
    let buf = BufWriter::new(file);
    let gz = path
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"));
    Ok(if gz {
        Box::new(GzEncoder::new(buf, Compression::default()))
    } else {
        Box::new(buf)
    })
}

fn write_record<W: Write>(w: &mut W, id: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
    w.write_all(b"@").rs_context("writing FASTQ")?;
    w.write_all(id).rs_context("writing FASTQ")?;
    w.write_all(b"\n").rs_context("writing FASTQ")?;
    w.write_all(seq).rs_context("writing FASTQ")?;
    w.write_all(b"\n+\n").rs_context("writing FASTQ")?;
    w.write_all(qual).rs_context("writing FASTQ")?;
    w.write_all(b"\n").rs_context("writing FASTQ")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_fixture(path: &Path, records: &[(&str, &str, &str)]) {
        use std::fmt::Write as _;
        let mut s = String::new();
        for (id, seq, qual) in records {
            writeln!(&mut s, "@{id}\n{seq}\n+\n{qual}").unwrap();
        }
        std::fs::write(path, s).unwrap();
    }

    #[test]
    fn full_seq_removes_exact_duplicates() {
        let tmp = tempfile::tempdir().unwrap();
        let inp = tmp.path().join("in.fq");
        let out = tmp.path().join("out.fq");
        write_fixture(
            &inp,
            &[
                ("r1", "ACGTACGT", "IIIIIIII"),
                ("r2", "ACGTACGT", "IIIIIIII"),
                ("r3", "TTTTAAAA", "IIIIIIII"),
            ],
        );
        let r = run_se(
            &inp,
            &out,
            DedupConfig {
                mode: DedupMode::FullSeq,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(r.reads_in, 3);
        assert_eq!(r.reads_out, 2);
        assert_eq!(r.duplicates_removed, 1);
        let content = std::fs::read_to_string(&out).unwrap();
        assert!(content.contains("@r1\n"));
        assert!(!content.contains("@r2\n"));
        assert!(content.contains("@r3\n"));
    }

    #[test]
    fn kmer_bin_groups_tail_matches() {
        let tmp = tempfile::tempdir().unwrap();
        let inp = tmp.path().join("in.fq");
        let out = tmp.path().join("out.fq");
        let common_tail = "AAACCCGGGTTT";
        write_fixture(
            &inp,
            &[
                ("r1", &format!("AAAAA{common_tail}"), &"I".repeat(17)),
                ("r2", &format!("CCCCC{common_tail}"), &"I".repeat(17)),
                ("r3", "GGGGGGGGGGGGGGGGG", &"I".repeat(17)),
            ],
        );
        let r = run_se(
            &inp,
            &out,
            DedupConfig {
                mode: DedupMode::KmerBin,
                k: 12,
                tail_offset: 0,
            },
        )
        .unwrap();
        assert_eq!(r.reads_in, 3);
        assert_eq!(r.reads_out, 2, "r2 should bin with r1");
        assert_eq!(r.duplicates_removed, 1);
    }

    #[test]
    fn short_read_falls_back_to_full_seq_hash() {
        let tmp = tempfile::tempdir().unwrap();
        let inp = tmp.path().join("in.fq");
        let out = tmp.path().join("out.fq");
        write_fixture(
            &inp,
            &[
                ("r1", "ACGT", "IIII"),
                ("r2", "ACGT", "IIII"),
                ("r3", "TTTT", "IIII"),
            ],
        );
        let r = run_se(
            &inp,
            &out,
            DedupConfig {
                mode: DedupMode::KmerBin,
                k: 12,
                tail_offset: 0,
            },
        )
        .unwrap();
        assert_eq!(r.reads_out, 2);
    }
}
