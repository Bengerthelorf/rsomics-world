/// Per-pass read statistics, accumulated on the fly.
///
/// `pre_filter` and `post_filter` are independent `ReadStats` instances so the
/// JSON report can show both `summary.before_filtering` and
/// `summary.after_filtering` views — matching fastp's report layout.
#[derive(Debug, Clone, Default)]
pub struct ReadStats {
    pub total_reads: u64,
    pub total_bases: u64,
    pub q20_bases: u64,
    pub q30_bases: u64,
    pub gc_bases: u64,
    pub n_bases: u64,
    pub cycles: Vec<CycleStat>,
}

/// Per-cycle counts. Index = 0-based position along the read. `qual_sum` is
/// the sum of Phred values at that cycle across all reads that reached it;
/// divide by the per-cycle total to get the mean quality.
#[derive(Debug, Clone, Default)]
pub struct CycleStat {
    pub count_a: u64,
    pub count_c: u64,
    pub count_g: u64,
    pub count_t: u64,
    pub count_n: u64,
    pub qual_sum: u64,
}

impl CycleStat {
    #[must_use]
    pub fn total(&self) -> u64 {
        self.count_a + self.count_c + self.count_g + self.count_t + self.count_n
    }
}

impl ReadStats {
    /// Fold a single read into the running totals.
    pub fn observe(&mut self, seq: &[u8], qual: &[u8]) {
        self.total_reads += 1;
        self.total_bases += seq.len() as u64;

        if self.cycles.len() < seq.len() {
            self.cycles.resize(seq.len(), CycleStat::default());
        }

        for (i, &b) in seq.iter().enumerate() {
            let c = &mut self.cycles[i];
            match b {
                b'A' | b'a' => c.count_a += 1,
                b'C' | b'c' => {
                    c.count_c += 1;
                    self.gc_bases += 1;
                }
                b'G' | b'g' => {
                    c.count_g += 1;
                    self.gc_bases += 1;
                }
                b'T' | b't' => c.count_t += 1,
                b'N' | b'n' => {
                    c.count_n += 1;
                    self.n_bases += 1;
                }
                _ => {}
            }
        }

        for (i, &q) in qual.iter().enumerate() {
            // Sanger encoding: q = Phred + 33.
            let phred = q.saturating_sub(33);
            self.cycles[i].qual_sum += u64::from(phred);
            if phred >= 20 {
                self.q20_bases += 1;
                if phred >= 30 {
                    self.q30_bases += 1;
                }
            }
        }
    }

    // f64 from u64 base counts only loses precision past 2^53 bases — that's
    // about 9 petabytes of FASTQ. Not a concern at any realistic input size.
    #[allow(clippy::cast_precision_loss)]
    #[must_use]
    pub fn q20_rate(&self) -> f64 {
        if self.total_bases == 0 {
            0.0
        } else {
            self.q20_bases as f64 / self.total_bases as f64
        }
    }

    #[allow(clippy::cast_precision_loss)]
    #[must_use]
    pub fn q30_rate(&self) -> f64 {
        if self.total_bases == 0 {
            0.0
        } else {
            self.q30_bases as f64 / self.total_bases as f64
        }
    }

    #[allow(clippy::cast_precision_loss)]
    #[must_use]
    pub fn gc_content(&self) -> f64 {
        if self.total_bases == 0 {
            0.0
        } else {
            self.gc_bases as f64 / self.total_bases as f64
        }
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    #[test]
    fn empty_stats_have_zero_rates() {
        let s = ReadStats::default();
        assert_eq!(s.total_reads, 0);
        assert_eq!(s.q20_rate(), 0.0);
        assert_eq!(s.q30_rate(), 0.0);
        assert_eq!(s.gc_content(), 0.0);
    }

    #[test]
    fn single_read_aggregates_correctly() {
        let mut s = ReadStats::default();
        // 4 bases, all G, all Q40 (I = byte 73, phred = 40).
        s.observe(b"GGGG", b"IIII");
        assert_eq!(s.total_reads, 1);
        assert_eq!(s.total_bases, 4);
        assert_eq!(s.gc_bases, 4);
        assert_eq!(s.n_bases, 0);
        assert_eq!(s.q20_bases, 4);
        assert_eq!(s.q30_bases, 4);
        assert!((s.gc_content() - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn mixed_quality_partitions_q20_q30() {
        let mut s = ReadStats::default();
        // 'I' = Q40, '5' = Q20, '!' = Q0.
        s.observe(b"ACGT", b"I5!I");
        assert_eq!(s.q20_bases, 3);
        assert_eq!(s.q30_bases, 2);
    }

    #[test]
    fn n_bases_excluded_from_gc() {
        let mut s = ReadStats::default();
        s.observe(b"NCGN", b"IIII");
        assert_eq!(s.n_bases, 2);
        assert_eq!(s.gc_bases, 2);
    }
}
