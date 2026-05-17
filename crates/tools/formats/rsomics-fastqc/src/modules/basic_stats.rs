use std::fmt::Write as _;

use super::{ModuleStatus, QcModule, Record};

/// `FastQC` "Basic Statistics" — never warns or fails (always Pass). The
/// `%GC` it reports is an integer over all called bases; encoding is
/// inferred from the lowest quality byte (clean-room: modern Phred+33 data
/// reports "Sanger / Illumina 1.9").
pub struct BasicStats {
    filename: String,
    n_seqs: u64,
    min_len: u64,
    max_len: u64,
    gc: u64,
    total_bases: u64,
    min_qual_byte: u8,
}

impl BasicStats {
    #[must_use]
    pub fn new(filename: &str) -> Self {
        Self {
            filename: filename.to_owned(),
            n_seqs: 0,
            min_len: u64::MAX,
            max_len: 0,
            gc: 0,
            total_bases: 0,
            min_qual_byte: u8::MAX,
        }
    }

    fn encoding(&self) -> &'static str {
        // FastQC reports Phred+33 data as "Sanger / Illumina 1.9"; only
        // legacy Phred+64 runs (lowest byte ≥ 64) read as "Illumina 1.3".
        if self.min_qual_byte < 64 {
            "Sanger / Illumina 1.9"
        } else {
            "Illumina 1.3"
        }
    }

    fn seq_len_field(&self) -> String {
        if self.min_len == self.max_len {
            self.max_len.to_string()
        } else {
            format!("{}-{}", self.min_len, self.max_len)
        }
    }
}

impl QcModule for BasicStats {
    fn name(&self) -> &'static str {
        "Basic Statistics"
    }

    fn process(&mut self, rec: &Record) {
        self.n_seqs += 1;
        let len = rec.seq.len() as u64;
        self.min_len = self.min_len.min(len);
        self.max_len = self.max_len.max(len);
        self.total_bases += len;
        for &b in rec.seq {
            if matches!(b, b'G' | b'C' | b'g' | b'c') {
                self.gc += 1;
            }
        }
        for &q in rec.qual {
            self.min_qual_byte = self.min_qual_byte.min(q);
        }
    }

    fn finalize(&mut self) {
        if self.n_seqs == 0 {
            self.min_len = 0;
        }
    }

    fn status(&self) -> ModuleStatus {
        ModuleStatus::Pass
    }

    fn write_data(&self, out: &mut String) {
        #[allow(clippy::cast_precision_loss)]
        let gc_pct = if self.total_bases == 0 {
            0
        } else {
            // 0..=100 percentage after round(), fits u64
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let v = (self.gc as f64 * 100.0 / self.total_bases as f64).round() as u64;
            v
        };
        out.push_str("#Measure\tValue\n");
        let _ = writeln!(out, "Filename\t{}", self.filename);
        out.push_str("File type\tConventional base calls\n");
        let _ = writeln!(out, "Encoding\t{}", self.encoding());
        let _ = writeln!(out, "Total Sequences\t{}", self.n_seqs);
        out.push_str("Sequences flagged as poor quality\t0\n");
        let _ = writeln!(out, "Sequence length\t{}", self.seq_len_field());
        let _ = writeln!(out, "%GC\t{gc_pct}");
    }
}
