use std::io::{self, Read};
use std::path::Path;

#[cfg(target_os = "linux")]
use std::fs::File;
#[cfg(target_os = "linux")]
use std::mem::MaybeUninit;

#[cfg(target_os = "linux")]
use isal_sys::igzip_lib::{
    ISAL_DECOMP_OK, ISAL_END_INPUT, ISAL_GZIP, ISAL_OUT_OVERFLOW, inflate_state,
    isal_block_state_ISAL_BLOCK_FINISH as ISAL_BLOCK_FINISH, isal_gzip_header,
    isal_gzip_header_init, isal_inflate, isal_inflate_init, isal_read_gzip_header,
};

#[cfg(target_os = "linux")]
const IGZIP_IN: usize = 1 << 22;
#[cfg(target_os = "linux")]
const FQ_BUF: usize = 1 << 23;

#[cfg(not(target_os = "linux"))]
pub struct GzReader;

#[cfg(not(target_os = "linux"))]
impl GzReader {
    pub fn new(_path: &Path) -> io::Result<Self> {
        Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "rsomics-igzip: the ISA-L igzip backend is Linux-only on this build",
        ))
    }
}

#[cfg(not(target_os = "linux"))]
impl Read for GzReader {
    fn read(&mut self, _buf: &mut [u8]) -> io::Result<usize> {
        unreachable!("rsomics-igzip GzReader cannot be constructed off Linux")
    }
}

#[cfg(target_os = "linux")]
pub struct GzReader {
    file: File,
    in_buf: Vec<u8>,
    out_buf: Vec<u8>,
    out_pos: usize,
    out_end: usize,
    // ~85 KiB — heap-boxed to avoid stack overflow.
    state: Box<inflate_state>,
    done: bool,
}

#[cfg(target_os = "linux")]
impl GzReader {
    pub fn new(path: &Path) -> io::Result<Self> {
        let mut file = File::open(path)
            .map_err(|e| io::Error::new(e.kind(), format!("igzip open {}: {e}", path.display())))?;

        // SAFETY: inflate_state is #[repr(C)] with only integer/pointer/byte-array
        // fields — the all-zero bit pattern is valid (no enum/niche/NonNull).
        // isal_inflate_init sets all scalar fields before ISA-L reads them.
        let mut state = unsafe {
            let mut s: Box<MaybeUninit<inflate_state>> = Box::new(MaybeUninit::zeroed());
            isal_inflate_init(s.as_mut_ptr().cast::<inflate_state>());
            s.assume_init()
        };

        state.crc_flag = ISAL_GZIP;

        let mut in_buf = vec![0u8; IGZIP_IN];
        let n = file.read(&mut in_buf)?;
        if n == 0 {
            return Ok(Self {
                file,
                in_buf,
                out_buf: vec![0u8; FQ_BUF],
                out_pos: 0,
                out_end: 0,
                state,
                done: true,
            });
        }

        // SAFETY: in_buf outlives the call; state is initialised; n ≤ in_buf.len().
        let ret = unsafe {
            let mut gz_hdr: MaybeUninit<isal_gzip_header> = MaybeUninit::uninit();
            isal_gzip_header_init(gz_hdr.as_mut_ptr());
            state.next_in = in_buf.as_mut_ptr();
            state.avail_in = n as u32;
            isal_read_gzip_header(&mut *state as *mut inflate_state, gz_hdr.as_mut_ptr())
        };

        if ret != ISAL_DECOMP_OK as i32 && ret != ISAL_END_INPUT as i32 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("igzip: gzip header parse failed (isal error {ret})"),
            ));
        }

        Ok(Self {
            file,
            in_buf,
            out_buf: vec![0u8; FQ_BUF],
            out_pos: 0,
            out_end: 0,
            state,
            done: false,
        })
    }

    fn refill(&mut self) -> io::Result<bool> {
        loop {
            if self.state.avail_in == 0 {
                let n = self.file.read(&mut self.in_buf)?;
                if n == 0 {
                    // ISAL_BLOCK_FINISH is set only after the gzip trailer checksum validates — the sole clean-EOF signal; any other state at n==0 means truncation
                    if self.state.block_state == ISAL_BLOCK_FINISH {
                        return Ok(false);
                    }
                    return Err(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "igzip: unexpected end of gzip stream (truncated or corrupt)",
                    ));
                }
                self.state.next_in = self.in_buf.as_mut_ptr();
                self.state.avail_in = n as u32;
            }

            self.state.next_out = self.out_buf.as_mut_ptr();
            self.state.avail_out = FQ_BUF as u32;

            // SAFETY: next_in/next_out point into self.in_buf/out_buf; avail_in/avail_out are set above.
            let ret = unsafe { isal_inflate(&mut *self.state as *mut inflate_state) };

            let produced = FQ_BUF - self.state.avail_out as usize;

            if produced > 0 {
                self.out_pos = 0;
                self.out_end = produced;
            }

            match ret as u32 {
                ISAL_DECOMP_OK => {
                    if produced > 0 {
                        return Ok(true);
                    }
                }
                ISAL_OUT_OVERFLOW => {
                    if produced > 0 {
                        return Ok(true);
                    }
                    return Err(io::Error::other(
                        "igzip: ISAL_OUT_OVERFLOW with zero output — decompressor bug",
                    ));
                }
                ISAL_END_INPUT => {
                    if produced > 0 {
                        return Ok(true);
                    }
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("igzip: decompression error (isal error {ret})"),
                    ));
                }
            }

            if self.state.block_state == ISAL_BLOCK_FINISH && self.state.avail_in > 0 {
                let saved_next_in = self.state.next_in;
                let saved_avail_in = self.state.avail_in;

                // SAFETY: state is valid; isal_inflate_init zero-resets all fields.
                unsafe {
                    isal_inflate_init(&mut *self.state as *mut inflate_state);
                }
                self.state.crc_flag = ISAL_GZIP;

                // SAFETY: saved_next_in points into in_buf (still alive); saved_avail_in ≤ in_buf.len().
                let hdr_ret = unsafe {
                    let mut gz_hdr: MaybeUninit<isal_gzip_header> = MaybeUninit::uninit();
                    isal_gzip_header_init(gz_hdr.as_mut_ptr());
                    self.state.next_in = saved_next_in;
                    self.state.avail_in = saved_avail_in;
                    isal_read_gzip_header(
                        &mut *self.state as *mut inflate_state,
                        gz_hdr.as_mut_ptr(),
                    )
                };

                if hdr_ret != ISAL_DECOMP_OK as i32 && hdr_ret != ISAL_END_INPUT as i32 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "igzip: concatenated member header parse failed (isal error {hdr_ret})"
                        ),
                    ));
                }

                if produced > 0 {
                    return Ok(true);
                }
            } else if self.state.block_state == ISAL_BLOCK_FINISH && self.state.avail_in == 0 {
                if produced > 0 {
                    return Ok(true);
                }
                return Ok(false);
            }
        }
    }
}

#[cfg(target_os = "linux")]
impl Read for GzReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.done {
            return Ok(0);
        }

        if self.out_pos < self.out_end {
            let available = self.out_end - self.out_pos;
            let to_copy = available.min(buf.len());
            buf[..to_copy].copy_from_slice(&self.out_buf[self.out_pos..self.out_pos + to_copy]);
            self.out_pos += to_copy;
            return Ok(to_copy);
        }

        match self.refill() {
            Ok(false) => {
                self.done = true;
                Ok(0)
            }
            Ok(true) => {
                let available = self.out_end - self.out_pos;
                let to_copy = available.min(buf.len());
                buf[..to_copy].copy_from_slice(&self.out_buf[self.out_pos..self.out_pos + to_copy]);
                self.out_pos += to_copy;
                Ok(to_copy)
            }
            Err(e) => {
                self.done = true;
                Err(e)
            }
        }
    }
}

#[cfg(all(test, target_os = "linux"))]
mod tests {
    use std::io::{Read, Write};

    use flate2::Compression;
    use flate2::write::GzEncoder;

    use super::*;

    fn gz_encode(data: &[u8]) -> Vec<u8> {
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(data).unwrap();
        enc.finish().unwrap()
    }

    fn gz_encode_multi(members: &[&[u8]]) -> Vec<u8> {
        let mut out = Vec::new();
        for m in members {
            out.extend(gz_encode(m));
        }
        out
    }

    fn flate2_decompress(gz: &[u8]) -> Vec<u8> {
        let mut dec = flate2::read::MultiGzDecoder::new(gz);
        let mut out = Vec::new();
        dec.read_to_end(&mut out).unwrap();
        out
    }

    fn write_tmp(bytes: &[u8]) -> tempfile::NamedTempFile {
        let mut f = tempfile::Builder::new().suffix(".gz").tempfile().unwrap();
        f.write_all(bytes).unwrap();
        f.flush().unwrap();
        f
    }

    fn igzip_decompress(path: &Path) -> io::Result<Vec<u8>> {
        let mut r = GzReader::new(path)?;
        let mut out = Vec::new();
        r.read_to_end(&mut out)?;
        Ok(out)
    }

    #[test]
    fn single_member_small() {
        let data = b"hello igzip world\n";
        let gz = gz_encode(data);
        let f = write_tmp(&gz);
        let got = igzip_decompress(f.path()).unwrap();
        assert_eq!(got, data, "single-member round-trip mismatch");
    }

    #[test]
    fn single_member_identical_to_flate2() {
        let data = b"@read1 desc\nACGTACGT\n+\nIIIIIIII\n@read2\nTTTT\n+\nFFFF\n";
        let gz = gz_encode(data);
        let f = write_tmp(&gz);
        let got = igzip_decompress(f.path()).unwrap();
        let expected = flate2_decompress(&gz);
        assert_eq!(got, expected, "igzip vs flate2 output mismatch");
    }

    #[test]
    fn multi_member_two_records() {
        let a = b"@r1\nACGT\n+\nIIII\n";
        let b = b"@r2\nTTTT\n+\nFFFF\n";
        let gz = gz_encode_multi(&[a, b]);
        let f = write_tmp(&gz);
        let got = igzip_decompress(f.path()).unwrap();
        let expected = flate2_decompress(&gz);
        assert_eq!(got, expected, "multi-member round-trip mismatch");
    }

    #[test]
    fn multi_member_identical_to_flate2() {
        let parts: Vec<Vec<u8>> = (0..5)
            .map(|i| format!("@read{i}\nACGTACGT\n+\nIIIIIIII\n").into_bytes())
            .collect();
        let slices: Vec<&[u8]> = parts.iter().map(|v| v.as_slice()).collect();
        let gz = gz_encode_multi(&slices);
        let f = write_tmp(&gz);
        let got = igzip_decompress(f.path()).unwrap();
        let expected = flate2_decompress(&gz);
        assert_eq!(got, expected);
    }

    #[test]
    fn data_larger_than_fq_buf() {
        let data: Vec<u8> = (0u8..=255).cycle().take(12 * 1024 * 1024).collect();
        let gz = gz_encode(&data);
        let f = write_tmp(&gz);
        let got = igzip_decompress(f.path()).unwrap();
        assert_eq!(got.len(), data.len(), "length mismatch for >FQ_BUF input");
        assert_eq!(got, data, "content mismatch for >FQ_BUF input");
    }

    #[test]
    fn empty_gzip_yields_empty_output() {
        let gz = gz_encode(b"");
        let f = write_tmp(&gz);
        let got = igzip_decompress(f.path()).unwrap();
        assert!(got.is_empty(), "empty gz must decompress to empty");
    }

    #[test]
    fn zero_byte_file_yields_empty_output() {
        let f = write_tmp(b"");
        let got = igzip_decompress(f.path()).unwrap();
        assert!(got.is_empty(), "zero-byte file must yield empty output");
    }

    #[test]
    fn truncated_stream_errors_loudly() {
        let data = b"@r1\nACGT\n+\nIIII\n";
        let mut gz = gz_encode(data);
        let new_len = gz.len().saturating_sub(8);
        gz.truncate(new_len);
        let f = write_tmp(&gz);
        let result = igzip_decompress(f.path());
        assert!(result.is_err(), "truncated gz must return Err, got Ok");
    }

    #[test]
    fn trailer_truncation_one_through_eight_bytes_all_error() {
        let data = b"@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n@r2\nTTTT\n+\nFFFF\n";
        let full = gz_encode(data);
        for k in 1..=8 {
            let mut gz = full.clone();
            gz.truncate(gz.len() - k);
            let f = write_tmp(&gz);
            assert!(
                igzip_decompress(f.path()).is_err(),
                "truncating {k} trailer byte(s) must Err, not silently truncate"
            );
        }
    }

    #[test]
    fn multi_member_large_spanning_input_buffer() {
        // xorshift64 pseudo-random ACGT compresses poorly so the stream really exceeds the 4 MiB input buffer
        let mut s: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut rng_acgt = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|_| {
                    s ^= s << 13;
                    s ^= s >> 7;
                    s ^= s << 17;
                    b"ACGT"[(s & 3) as usize]
                })
                .collect()
        };
        let members: Vec<Vec<u8>> = (0..50)
            .map(|_| {
                let mut m = Vec::new();
                for i in 0..3000usize {
                    let seq = rng_acgt(150);
                    m.extend_from_slice(format!("@read{i}\n").as_bytes());
                    m.extend_from_slice(&seq);
                    m.extend_from_slice(b"\n+\n");
                    m.extend_from_slice(&[b'I'; 150]);
                    m.push(b'\n');
                }
                m
            })
            .collect();
        let slices: Vec<&[u8]> = members.iter().map(Vec::as_slice).collect();
        let gz = gz_encode_multi(&slices);
        assert!(
            gz.len() > 4 * (1 << 20),
            "fixture must exceed the 4 MiB input buffer (got {})",
            gz.len()
        );
        let f = write_tmp(&gz);
        let got = igzip_decompress(f.path()).unwrap();
        let expected = flate2_decompress(&gz);
        assert_eq!(
            got.len(),
            expected.len(),
            "length mismatch — member dropped?"
        );
        assert_eq!(got, expected, "multi-member spanning-input-buffer mismatch");
    }

    #[test]
    fn corrupt_header_errors_loudly() {
        let mut gz = gz_encode(b"some data");
        gz[0] = 0xFF;
        gz[1] = 0xFF;
        let f = write_tmp(&gz);
        let result = GzReader::new(f.path());
        assert!(
            result.is_err(),
            "corrupt gzip header must error at construction"
        );
    }

    #[test]
    fn fastq_fixture_round_trip() {
        let fq: Vec<u8> = (0..100)
            .flat_map(|i| format!("@read{i} desc{i}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n").into_bytes())
            .collect();
        let gz = gz_encode(&fq);
        let f = write_tmp(&gz);
        let got = igzip_decompress(f.path()).unwrap();
        assert_eq!(got, fq, "FASTQ fixture round-trip failed");
    }
}
