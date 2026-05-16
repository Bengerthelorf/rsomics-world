#![allow(clippy::missing_errors_doc, clippy::missing_panics_doc)]

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Genotype {
    /// Homozygous for allele 1 (PLINK code 0b00).
    HomA1,
    /// Missing genotype (PLINK code 0b01).
    Missing,
    /// Heterozygous (PLINK code 0b10).
    Het,
    /// Homozygous for allele 2 (PLINK code 0b11).
    HomA2,
}

#[derive(Debug, Clone)]
pub struct Variant {
    pub chrom: String,
    pub id: String,
    pub cm: f64,
    pub pos: u64,
    pub a1: String,
    pub a2: String,
}

#[derive(Debug, Clone)]
pub struct Sample {
    pub fid: String,
    pub iid: String,
    pub pid: String,
    pub mid: String,
    pub sex: u8,
    pub phen: String,
}

#[derive(Debug, Clone)]
pub struct Pgen {
    pub variants: Vec<Variant>,
    pub samples: Vec<Sample>,
    /// Variant-major matrix: rows = variants, cols = samples. Length =
    /// `n_variants × n_samples`; access via `gt[v_idx * n_samples + s_idx]`.
    pub gt: Vec<Genotype>,
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum PgenError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("bad .bed magic: got [{0:#x}, {1:#x}, {2:#x}], expected [0x6c, 0x1b, 0x01]")]
    BadMagic(u8, u8, u8),
    #[error(
        "sample-major .bed (3rd byte = 0x00) — not supported in 0.1; rerun PLINK with --make-bed"
    )]
    SampleMajor,
    #[error("malformed {file} line {line}: {reason}")]
    Malformed {
        file: String,
        line: usize,
        reason: String,
    },
    #[error(
        ".bed size mismatch: have {have} bytes, expected {expected} for {n_variants}×{n_samples}"
    )]
    SizeMismatch {
        have: u64,
        expected: u64,
        n_variants: usize,
        n_samples: usize,
    },
}

pub type Result<T> = std::result::Result<T, PgenError>;

impl Pgen {
    /// Load a PLINK1 fileset by `prefix` (`<prefix>.bim`, `.fam`, `.bed`).
    pub fn load(prefix: &Path) -> Result<Self> {
        let bim = prefix.with_extension("bim");
        let fam = prefix.with_extension("fam");
        let bed = prefix.with_extension("bed");
        let variants = parse_bim(&bim)?;
        let samples = parse_fam(&fam)?;
        let gt = parse_bed(&bed, variants.len(), samples.len())?;
        Ok(Self {
            variants,
            samples,
            gt,
        })
    }

    #[must_use]
    pub fn n_variants(&self) -> usize {
        self.variants.len()
    }
    #[must_use]
    pub fn n_samples(&self) -> usize {
        self.samples.len()
    }

    /// Genotype call for `(variant_idx, sample_idx)`.
    ///
    /// # Panics
    /// Panics on out-of-range indices.
    #[must_use]
    pub fn get(&self, v: usize, s: usize) -> Genotype {
        self.gt[v * self.samples.len() + s]
    }
}

fn parse_bim(path: &Path) -> Result<Vec<Variant>> {
    let f = File::open(path)?;
    let mut out = Vec::new();
    for (lineno, line) in BufReader::new(f).lines().enumerate() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 6 {
            return Err(PgenError::Malformed {
                file: path.display().to_string(),
                line: lineno + 1,
                reason: format!("{} fields, expected 6", fields.len()),
            });
        }
        let cm: f64 = fields[2].parse().map_err(|_| PgenError::Malformed {
            file: path.display().to_string(),
            line: lineno + 1,
            reason: format!("bad cM {:?}", fields[2]),
        })?;
        let pos: u64 = fields[3].parse().map_err(|_| PgenError::Malformed {
            file: path.display().to_string(),
            line: lineno + 1,
            reason: format!("bad pos {:?}", fields[3]),
        })?;
        out.push(Variant {
            chrom: fields[0].to_string(),
            id: fields[1].to_string(),
            cm,
            pos,
            a1: fields[4].to_string(),
            a2: fields[5].to_string(),
        });
    }
    Ok(out)
}

fn parse_fam(path: &Path) -> Result<Vec<Sample>> {
    let f = File::open(path)?;
    let mut out = Vec::new();
    for (lineno, line) in BufReader::new(f).lines().enumerate() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 6 {
            return Err(PgenError::Malformed {
                file: path.display().to_string(),
                line: lineno + 1,
                reason: format!("{} fields, expected 6", fields.len()),
            });
        }
        let sex: u8 = fields[4].parse().map_err(|_| PgenError::Malformed {
            file: path.display().to_string(),
            line: lineno + 1,
            reason: format!("bad sex code {:?}", fields[4]),
        })?;
        out.push(Sample {
            fid: fields[0].to_string(),
            iid: fields[1].to_string(),
            pid: fields[2].to_string(),
            mid: fields[3].to_string(),
            sex,
            phen: fields[5].to_string(),
        });
    }
    Ok(out)
}

fn parse_bed(path: &Path, n_variants: usize, n_samples: usize) -> Result<Vec<Genotype>> {
    let mut f = File::open(path)?;
    let mut header = [0_u8; 3];
    f.read_exact(&mut header)?;
    if header[0] != 0x6c || header[1] != 0x1b {
        return Err(PgenError::BadMagic(header[0], header[1], header[2]));
    }
    match header[2] {
        0x01 => {} // variant-major (only supported case)
        0x00 => return Err(PgenError::SampleMajor),
        other => return Err(PgenError::BadMagic(header[0], header[1], other)),
    }
    let bytes_per_variant = n_samples.div_ceil(4);
    let expected = (bytes_per_variant * n_variants) as u64;
    let have = f.metadata()?.len() - 3;
    if have != expected {
        return Err(PgenError::SizeMismatch {
            have,
            expected,
            n_variants,
            n_samples,
        });
    }
    let mut buf = vec![0_u8; bytes_per_variant];
    let mut out = Vec::with_capacity(n_variants * n_samples);
    for _ in 0..n_variants {
        f.read_exact(&mut buf)?;
        for s in 0..n_samples {
            let byte = buf[s / 4];
            let code = (byte >> ((s % 4) * 2)) & 0b11; // PLINK 2-bit little-endian
            out.push(match code {
                0b00 => Genotype::HomA1,
                0b01 => Genotype::Missing,
                0b10 => Genotype::Het,
                0b11 => Genotype::HomA2,
                _ => unreachable!(),
            });
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn tiny_fixture() -> (tempfile::TempDir, std::path::PathBuf) {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("toy");
        let mut bim = File::create(prefix.with_extension("bim")).unwrap();
        writeln!(bim, "1\trs1\t0\t100\tA\tG").unwrap();
        writeln!(bim, "1\trs2\t0\t200\tC\tT").unwrap();
        writeln!(bim, "2\trs3\t0\t300\tG\tC").unwrap();
        let mut fam = File::create(prefix.with_extension("fam")).unwrap();
        writeln!(fam, "F1\tI1\t0\t0\t1\t-9").unwrap();
        writeln!(fam, "F2\tI2\t0\t0\t2\t-9").unwrap();
        writeln!(fam, "F3\tI3\t0\t0\t1\t-9").unwrap();
        writeln!(fam, "F4\tI4\t0\t0\t2\t-9").unwrap();
        // .bed magic + 3 bytes (ceil(4 samples / 4) = 1 byte per variant).
        // Variant 0: [HomA1=00, Het=10, HomA2=11, Missing=01] packed
        //   little-endian 2-bit: byte = 0b01_11_10_00 = 0x78.
        // Variant 1: all HomA2 (0b11) → 0xFF.
        // Variant 2: all HomA1 (0b00) → 0x00.
        let mut bed = File::create(prefix.with_extension("bed")).unwrap();
        bed.write_all(&[0x6c, 0x1b, 0x01, 0x78, 0xFF, 0x00])
            .unwrap();
        (dir, prefix)
    }

    #[test]
    fn load_tiny_fixture() {
        let (_d, prefix) = tiny_fixture();
        let pgen = Pgen::load(&prefix).unwrap();
        assert_eq!(pgen.n_variants(), 3);
        assert_eq!(pgen.n_samples(), 4);
        assert_eq!(pgen.get(0, 0), Genotype::HomA1);
        assert_eq!(pgen.get(0, 1), Genotype::Het);
        assert_eq!(pgen.get(0, 2), Genotype::HomA2);
        assert_eq!(pgen.get(0, 3), Genotype::Missing);
        for s in 0..4 {
            assert_eq!(pgen.get(1, s), Genotype::HomA2);
        }
        for s in 0..4 {
            assert_eq!(pgen.get(2, s), Genotype::HomA1);
        }
    }

    #[test]
    fn bad_magic_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("toy");
        let _ = File::create(prefix.with_extension("bim")).unwrap();
        let _ = File::create(prefix.with_extension("fam")).unwrap();
        let mut bed = File::create(prefix.with_extension("bed")).unwrap();
        bed.write_all(&[0x00, 0x00, 0x00]).unwrap();
        assert!(matches!(Pgen::load(&prefix), Err(PgenError::BadMagic(..))));
    }

    #[test]
    fn sample_major_rejected_with_clear_error() {
        let (_d, prefix) = tiny_fixture();
        let mut bed = File::create(prefix.with_extension("bed")).unwrap();
        bed.write_all(&[0x6c, 0x1b, 0x00, 0x78, 0xFF, 0x00])
            .unwrap();
        assert!(matches!(Pgen::load(&prefix), Err(PgenError::SampleMajor)));
    }

    #[test]
    fn bim_field_count_validated() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("toy");
        let mut bim = File::create(prefix.with_extension("bim")).unwrap();
        writeln!(bim, "1\trs1\t100").unwrap();
        File::create(prefix.with_extension("fam")).unwrap();
        File::create(prefix.with_extension("bed")).unwrap();
        let err = Pgen::load(&prefix).unwrap_err();
        assert!(matches!(err, PgenError::Malformed { .. }));
    }

    #[test]
    fn size_mismatch_rejected() {
        let (_d, prefix) = tiny_fixture();
        let mut bed = File::create(prefix.with_extension("bed")).unwrap();
        bed.write_all(&[0x6c, 0x1b, 0x01, 0x78]).unwrap();
        assert!(matches!(
            Pgen::load(&prefix),
            Err(PgenError::SizeMismatch { .. })
        ));
    }
}
