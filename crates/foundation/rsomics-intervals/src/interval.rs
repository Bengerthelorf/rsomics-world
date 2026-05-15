use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    #[must_use]
    pub const fn as_byte(self) -> u8 {
        match self {
            Self::Forward => b'+',
            Self::Reverse => b'-',
        }
    }

    #[must_use]
    pub const fn from_byte(b: u8) -> Option<Self> {
        match b {
            b'+' => Some(Self::Forward),
            b'-' => Some(Self::Reverse),
            _ => None,
        }
    }
}

/// 0-based, half-open `[start, end)` interval on a named contig.
///
/// Intentionally tiny (chrom + range + optional strand) — hot-path algorithms
/// iterate over millions of these. Per-record extras (name, score, BED12
/// blocks) belong on a wrapper that owns an `Interval`.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Interval {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub strand: Option<Strand>,
}

impl Interval {
    #[allow(clippy::missing_errors_doc)]
    pub fn new(chrom: impl Into<String>, start: u64, end: u64) -> Result<Self, IntervalError> {
        if start >= end {
            return Err(IntervalError::Empty { start, end });
        }
        Ok(Self {
            chrom: chrom.into(),
            start,
            end,
            strand: None,
        })
    }

    #[allow(clippy::missing_errors_doc)]
    pub fn with_strand(
        chrom: impl Into<String>,
        start: u64,
        end: u64,
        strand: Strand,
    ) -> Result<Self, IntervalError> {
        let mut iv = Self::new(chrom, start, end)?;
        iv.strand = Some(strand);
        Ok(iv)
    }

    #[must_use]
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    /// Always false — `Interval::new` rejects empty intervals at construction.
    /// Present so generic code can call `.is_empty()` without a special case.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        false
    }

    #[must_use]
    pub fn overlaps(&self, other: &Self) -> bool {
        self.chrom == other.chrom && self.start < other.end && other.start < self.end
    }

    #[must_use]
    pub fn contains(&self, other: &Self) -> bool {
        self.chrom == other.chrom && self.start <= other.start && other.end <= self.end
    }

    #[must_use]
    pub fn overlap_bases(&self, other: &Self) -> u64 {
        if self.chrom != other.chrom {
            return 0;
        }
        let lo = self.start.max(other.start);
        let hi = self.end.min(other.end);
        hi.saturating_sub(lo)
    }
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum IntervalError {
    #[error("empty or inverted interval: start={start} >= end={end}")]
    Empty { start: u64, end: u64 },
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(chrom: &str, start: u64, end: u64) -> Interval {
        Interval::new(chrom, start, end).unwrap()
    }

    #[test]
    fn empty_interval_rejected() {
        assert!(matches!(
            Interval::new("chr1", 100, 100),
            Err(IntervalError::Empty { .. })
        ));
        assert!(matches!(
            Interval::new("chr1", 200, 100),
            Err(IntervalError::Empty { .. })
        ));
    }

    #[test]
    fn length_is_end_minus_start() {
        assert_eq!(iv("chr1", 100, 150).len(), 50);
    }

    #[test]
    fn overlaps_same_chrom() {
        let a = iv("chr1", 100, 200);
        assert!(a.overlaps(&iv("chr1", 150, 250)));
        assert!(a.overlaps(&iv("chr1", 50, 150)));
        assert!(a.overlaps(&iv("chr1", 100, 200)));
        assert!(!a.overlaps(&iv("chr1", 200, 300)), "half-open touching");
        assert!(!a.overlaps(&iv("chr1", 0, 100)), "half-open touching low");
        assert!(!a.overlaps(&iv("chr1", 250, 300)));
    }

    #[test]
    fn overlaps_different_chrom_is_false() {
        assert!(!iv("chr1", 100, 200).overlaps(&iv("chr2", 100, 200)));
    }

    #[test]
    fn contains_is_inclusive_at_edges() {
        let outer = iv("chr1", 100, 200);
        assert!(outer.contains(&iv("chr1", 100, 200)));
        assert!(outer.contains(&iv("chr1", 120, 180)));
        assert!(!outer.contains(&iv("chr1", 100, 201)));
        assert!(!outer.contains(&iv("chr1", 99, 150)));
    }

    #[test]
    fn overlap_bases_counts_intersection_length() {
        assert_eq!(
            iv("chr1", 100, 200).overlap_bases(&iv("chr1", 150, 250)),
            50
        );
        assert_eq!(iv("chr1", 100, 200).overlap_bases(&iv("chr1", 200, 300)), 0);
        assert_eq!(iv("chr1", 100, 200).overlap_bases(&iv("chr2", 100, 200)), 0);
        assert_eq!(
            iv("chr1", 100, 200).overlap_bases(&iv("chr1", 120, 180)),
            60
        );
    }

    #[test]
    fn strand_round_trips() {
        assert_eq!(Strand::Forward.as_byte(), b'+');
        assert_eq!(Strand::from_byte(b'+'), Some(Strand::Forward));
        assert_eq!(Strand::from_byte(b'.'), None);
    }
}
