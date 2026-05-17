// Slab invariant: each slab the producer sends holds whole FASTQ records (4-line
// counting, immune to @/+ in quality). Carry holds the tail into the next block;
// a record larger than OUT_BUF accumulates until a boundary.
use std::collections::VecDeque;
use std::io::{BufReader, Cursor, Read};
use std::path::Path;
use std::thread;

use crossbeam_channel::{Receiver, bounded};
use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};

use crate::OwnedRecord;
use crate::parse::parse_record;

#[cfg(not(target_os = "linux"))]
const IN_BUF: usize = 4 * 1024 * 1024;
const OUT_BUF: usize = 8 * 1024 * 1024;
const CHAN_DEPTH: usize = 32;

pub struct GzReader {
    rx: Receiver<Result<Vec<u8>>>,
    parsed: VecDeque<OwnedRecord>,
    done: bool,
    pending_err: Option<RsomicsError>,
}

impl GzReader {
    pub fn open(path: &Path) -> Result<Self> {
        let path = path.to_path_buf();
        let (tx, rx) = bounded(CHAN_DEPTH);

        thread::Builder::new()
            .name(format!("seqio-gz-{}", path.display()))
            .spawn(move || producer(&path, &tx))
            .map_err(RsomicsError::Io)?;

        Ok(Self {
            rx,
            parsed: VecDeque::new(),
            done: false,
            pending_err: None,
        })
    }

    fn refill(&mut self) -> bool {
        loop {
            match self.rx.recv() {
                Err(_) => {
                    self.done = true;
                    return false;
                }
                Ok(Err(e)) => {
                    self.done = true;
                    self.pending_err = Some(e);
                    return false;
                }
                Ok(Ok(slab)) => {
                    if slab.is_empty() {
                        self.done = true;
                        return false;
                    }
                    match parse_slab_parallel(&slab) {
                        Err(e) => {
                            self.done = true;
                            self.pending_err = Some(e);
                            return false;
                        }
                        Ok(recs) => {
                            if recs.is_empty() {
                                continue;
                            }
                            self.parsed = recs.into();
                            return true;
                        }
                    }
                }
            }
        }
    }
}

impl Iterator for GzReader {
    type Item = Result<OwnedRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(e) = self.pending_err.take() {
            return Some(Err(e));
        }
        if self.done {
            return None;
        }
        if let Some(rec) = self.parsed.pop_front() {
            return Some(Ok(rec));
        }
        if self.refill() {
            self.parsed.pop_front().map(Ok)
        } else {
            self.pending_err.take().map(Err)
        }
    }
}

fn parse_slab_parallel(slab: &[u8]) -> Result<Vec<OwnedRecord>> {
    let starts = record_start_offsets(slab);
    if starts.is_empty() {
        return Ok(Vec::new());
    }

    let ranges: Vec<(usize, usize)> = starts
        .windows(2)
        .map(|w| (w[0], w[1]))
        .chain(std::iter::once((*starts.last().unwrap(), slab.len())))
        .collect();

    ranges
        .into_par_iter()
        .map(|(start, end)| {
            let slice = &slab[start..end];
            let mut cur = Cursor::new(slice);
            parse_record(&mut cur)?
                .ok_or_else(|| RsomicsError::InvalidInput("empty record slice in slab".into()))
        })
        .collect()
}

fn record_start_offsets(slab: &[u8]) -> Vec<usize> {
    if slab.is_empty() {
        return Vec::new();
    }
    let mut starts = vec![0usize];
    let mut newline_count = 0u8;
    for (i, &b) in slab.iter().enumerate() {
        if b == b'\n' {
            newline_count += 1;
            if newline_count == 4 {
                newline_count = 0;
                let next = i + 1;
                if next < slab.len() {
                    starts.push(next);
                }
            }
        }
    }
    starts
}

fn producer(path: &Path, tx: &crossbeam_channel::Sender<Result<Vec<u8>>>) {
    // without catch_unwind a panicking decode thread drops tx — consumer's recv() reads that as clean EOF, silently truncating
    let outcome =
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| produce_inner(path, tx)));
    match outcome {
        Ok(Ok(())) => {}
        Ok(Err(e)) => {
            let _ = tx.send(Err(e));
        }
        Err(_) => {
            let _ = tx.send(Err(RsomicsError::Io(std::io::Error::other(
                "seqio: gz decode thread panicked",
            ))));
        }
    }
}

fn produce_inner(path: &Path, tx: &crossbeam_channel::Sender<Result<Vec<u8>>>) -> Result<()> {
    let decoder = build_decoder(path)?;
    let mut rdr = BufReader::with_capacity(OUT_BUF, decoder);

    let mut carry: Vec<u8> = Vec::new();
    let mut block = vec![0u8; OUT_BUF];

    loop {
        let n = read_full_block(&mut rdr, &mut block)?;
        if n == 0 {
            if !carry.is_empty() && tx.send(Ok(carry)).is_err() {
                return Ok(());
            }
            let _ = tx.send(Ok(Vec::new()));
            return Ok(());
        }

        let mut candidate = carry;
        candidate.extend_from_slice(&block[..n]);

        match last_record_boundary(&candidate) {
            None => {
                carry = candidate;
            }
            Some(boundary) => {
                let slab = candidate[..boundary].to_vec();
                carry = candidate[boundary..].to_vec();

                if !slab.is_empty() && tx.send(Ok(slab)).is_err() {
                    return Ok(());
                }
            }
        }
    }
}

fn read_full_block<R: Read>(rdr: &mut R, buf: &mut [u8]) -> Result<usize> {
    let mut total = 0;
    while total < buf.len() {
        match rdr.read(&mut buf[total..]) {
            Ok(0) => break,
            Ok(n) => total += n,
            Err(e) => return Err(RsomicsError::Io(e)),
        }
    }
    Ok(total)
}

fn last_record_boundary(data: &[u8]) -> Option<usize> {
    let mut newline_count = 0u8;
    let mut last_boundary = None;
    for (i, &b) in data.iter().enumerate() {
        if b == b'\n' {
            newline_count += 1;
            if newline_count == 4 {
                newline_count = 0;
                last_boundary = Some(i + 1);
            }
        }
    }
    last_boundary
}

#[cfg(target_os = "linux")]
fn build_decoder(path: &Path) -> Result<Box<dyn Read>> {
    let reader = rsomics_igzip::GzReader::new(path).map_err(RsomicsError::Io)?;
    Ok(Box::new(reader))
}

#[cfg(not(target_os = "linux"))]
fn build_decoder(path: &Path) -> Result<Box<dyn Read>> {
    use std::fs::File;
    use std::io::BufReader as StdBuf;

    let file = File::open(path).map_err(|e| {
        RsomicsError::Io(std::io::Error::new(
            e.kind(),
            format!("opening {}: {e}", path.display()),
        ))
    })?;
    let raw = StdBuf::with_capacity(IN_BUF, file);
    Ok(Box::new(flate2::read::MultiGzDecoder::new(raw)))
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    fn make_plain_gz(content: &[u8]) -> Vec<u8> {
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(content).unwrap();
        enc.finish().unwrap()
    }

    fn write_tmp_gz(content: &[u8]) -> tempfile::NamedTempFile {
        let mut f = tempfile::Builder::new()
            .suffix(".fq.gz")
            .tempfile()
            .unwrap();
        f.write_all(&make_plain_gz(content)).unwrap();
        f.flush().unwrap();
        f
    }

    fn collect_gz(content: &[u8]) -> Vec<OwnedRecord> {
        let f = write_tmp_gz(content);
        GzReader::open(f.path())
            .unwrap()
            .collect::<Result<Vec<_>>>()
            .unwrap()
    }

    #[test]
    fn gz_round_trip_two_records() {
        let fq = b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nFFFF\n";
        let recs = collect_gz(fq);
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, b"r1");
        assert_eq!(recs[0].seq, b"ACGT");
        assert_eq!(recs[1].id, b"r2");
        assert_eq!(recs[1].seq, b"TTTT");
    }

    #[test]
    fn gz_empty_input_yields_no_records() {
        let recs = collect_gz(b"");
        assert!(recs.is_empty());
    }

    #[test]
    fn gz_large_batch_multi_slab() {
        let mut content = Vec::new();
        for i in 0..200_001usize {
            writeln!(
                content,
                "@r{i}\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII"
            )
            .unwrap();
        }
        let recs = collect_gz(&content);
        assert_eq!(recs.len(), 200_001);
    }

    #[test]
    fn gz_records_identical_to_plain_parse() {
        use std::io::Cursor;

        let fq: &[u8] = b"\
@read1 desc1\nACGTACGT\n+\nIIIIIIII\n\
@read2\nTTTTCCCC\n+\nFFFFAAAA\n\
@read3\nGGGGAAAA\n+\nHHHHHHHH\n";

        let expected: Vec<_> = {
            let mut cur = Cursor::new(fq);
            let mut out = Vec::new();
            while let Some(rec) = crate::parse::parse_record(&mut cur).unwrap() {
                out.push(rec);
            }
            out
        };

        let got = collect_gz(fq);
        assert_eq!(got.len(), expected.len(), "record count mismatch");
        for (i, (g, e)) in got.iter().zip(expected.iter()).enumerate() {
            assert_eq!(g.id, e.id, "id mismatch at record {i}");
            assert_eq!(g.seq, e.seq, "seq mismatch at record {i}");
            assert_eq!(g.qual, e.qual, "qual mismatch at record {i}");
        }
    }

    #[test]
    fn gz_truncated_data_errors_loudly() {
        let fq = b"@r1\nACGT\n+\nIIII\n";
        let mut gz = make_plain_gz(fq);
        let new_len = gz.len().saturating_sub(6);
        gz.truncate(new_len);
        let mut f = tempfile::Builder::new()
            .suffix(".fq.gz")
            .tempfile()
            .unwrap();
        f.write_all(&gz).unwrap();
        f.flush().unwrap();

        let result: Result<Vec<_>> = GzReader::open(f.path()).unwrap().collect();
        assert!(result.is_err(), "truncated gz must error loudly");
    }

    #[test]
    fn gz_missing_final_newline_last_record_intact() {
        let recs = collect_gz(b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nFFFF");
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[1].id, b"r2");
        assert_eq!(recs[1].seq, b"TTTT");
        assert_eq!(recs[1].qual, b"FFFF");
    }

    #[test]
    fn gz_record_split_across_members() {
        let m1 = make_plain_gz(b"@r1\nACGT\n+\nII");
        let m2 = make_plain_gz(b"II\n@r2\nTTTT\n+\nFFFF\n");
        let mut concatenated = m1;
        concatenated.extend_from_slice(&m2);
        let mut f = tempfile::Builder::new()
            .suffix(".fq.gz")
            .tempfile()
            .unwrap();
        f.write_all(&concatenated).unwrap();
        f.flush().unwrap();

        let recs: Vec<_> = GzReader::open(f.path())
            .unwrap()
            .collect::<Result<Vec<_>>>()
            .unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, b"r1");
        assert_eq!(recs[0].qual, b"IIII");
        assert_eq!(recs[1].id, b"r2");
        assert_eq!(recs[1].seq, b"TTTT");
    }

    #[test]
    fn slab_boundary_carry_two_records() {
        let long_seq: Vec<u8> = b"ACGT"
            .iter()
            .cycle()
            .copied()
            .take(3 * 1024 * 1024)
            .collect();
        let long_qual: Vec<u8> = b"I".iter().cycle().copied().take(3 * 1024 * 1024).collect();

        let mut fq = Vec::new();
        fq.extend_from_slice(b"@rec1\n");
        fq.extend_from_slice(&long_seq);
        fq.extend_from_slice(b"\n+\n");
        fq.extend_from_slice(&long_qual);
        fq.extend_from_slice(b"\n@rec2\n");
        fq.extend_from_slice(&long_seq);
        fq.extend_from_slice(b"\n+\n");
        fq.extend_from_slice(&long_qual);
        fq.extend_from_slice(b"\n");

        let recs = collect_gz(&fq);
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, b"rec1");
        assert_eq!(recs[0].seq, long_seq.as_slice());
        assert_eq!(recs[0].qual, long_qual.as_slice());
        assert_eq!(recs[1].id, b"rec2");
        assert_eq!(recs[1].seq, long_seq.as_slice());
        assert_eq!(recs[1].qual, long_qual.as_slice());
    }

    #[test]
    fn single_record_larger_than_out_buf() {
        let big_seq: Vec<u8> = b"ACGT"
            .iter()
            .cycle()
            .copied()
            .take(10 * 1024 * 1024)
            .collect();
        let big_qual: Vec<u8> = b"I"
            .iter()
            .cycle()
            .copied()
            .take(10 * 1024 * 1024)
            .collect();

        let mut fq = Vec::new();
        fq.extend_from_slice(b"@big\n");
        fq.extend_from_slice(&big_seq);
        fq.extend_from_slice(b"\n+\n");
        fq.extend_from_slice(&big_qual);
        fq.extend_from_slice(b"\n");

        let recs = collect_gz(&fq);
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].id, b"big");
        assert_eq!(recs[0].seq.len(), big_seq.len());
        assert_eq!(recs[0].qual.len(), big_qual.len());
    }

    #[test]
    fn crlf_records_through_slab() {
        let fq = b"@r1\r\nACGT\r\n+\r\nIIII\r\n@r2\r\nTTTT\r\n+\r\nFFFF\r\n";
        let recs = collect_gz(fq);
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, b"r1");
        assert_eq!(recs[0].seq, b"ACGT");
        assert_eq!(recs[1].id, b"r2");
        assert_eq!(recs[1].seq, b"TTTT");
    }

    #[test]
    fn at_in_quality_does_not_split_record() {
        let fq = b"@r1\nACGT\n+\n@III\n@r2\nTTTT\n+\nFFFF\n";
        let recs = collect_gz(fq);
        assert_eq!(
            recs.len(),
            2,
            "@ in quality must not create a spurious record"
        );
        assert_eq!(recs[0].id, b"r1");
        assert_eq!(recs[0].qual, b"@III");
        assert_eq!(recs[1].id, b"r2");
    }

    #[test]
    fn last_record_boundary_two_records() {
        let data = b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nFFFF\n";
        let pos = last_record_boundary(data).unwrap();
        assert_eq!(pos, data.len(), "boundary should be at end of both records");
    }

    #[test]
    fn last_record_boundary_one_complete_one_partial() {
        let data = b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nFF";
        let pos = last_record_boundary(data).unwrap();
        assert_eq!(&data[..pos], b"@r1\nACGT\n+\nIIII\n");
        assert_eq!(&data[pos..], b"@r2\nTTTT\n+\nFF");
    }

    #[test]
    fn last_record_boundary_no_complete_record() {
        let data = b"@r1\nACGT\n+\n";
        assert!(last_record_boundary(data).is_none());
    }

    #[test]
    fn last_record_boundary_empty() {
        assert!(last_record_boundary(b"").is_none());
    }

    #[test]
    fn record_start_offsets_two_records() {
        let slab = b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nFFFF\n";
        let offsets = record_start_offsets(slab);
        assert_eq!(offsets, vec![0, 16]);
        assert_eq!(slab[16], b'@');
    }

    #[test]
    fn record_start_offsets_single() {
        let slab = b"@r1\nACGT\n+\nIIII\n";
        assert_eq!(record_start_offsets(slab), vec![0]);
    }

    #[test]
    fn record_start_offsets_empty() {
        assert!(record_start_offsets(b"").is_empty());
    }

    #[test]
    fn record_start_offsets_at_in_quality() {
        let slab = b"@r1\nACGT\n+\n@III\n@r2\nTTTT\n+\nFFFF\n";
        let offsets = record_start_offsets(slab);
        assert_eq!(offsets, vec![0, 16]);
    }
}
