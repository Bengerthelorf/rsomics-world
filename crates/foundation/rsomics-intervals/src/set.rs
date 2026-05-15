
use std::collections::BTreeMap;

use crate::interval::Interval;

/// Chromosome-grouped collection of intervals. `BTreeMap` gives lexicographic
/// chrom iteration order — matches `bedtools sort` default.
#[derive(Debug, Default, Clone)]
pub struct IntervalSet {
    by_chrom: BTreeMap<String, Vec<Interval>>,
    sorted: bool,
}

impl IntervalSet {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    pub fn push(&mut self, interval: Interval) {
        self.by_chrom
            .entry(interval.chrom.clone())
            .or_default()
            .push(interval);
        self.sorted = false;
    }

    pub fn extend<I: IntoIterator<Item = Interval>>(&mut self, items: I) {
        for iv in items {
            self.by_chrom.entry(iv.chrom.clone()).or_default().push(iv);
        }
        self.sorted = false;
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.by_chrom.values().map(Vec::len).sum()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn sort(&mut self) {
        if self.sorted {
            return;
        }
        for ivs in self.by_chrom.values_mut() {
            ivs.sort_by(|a, b| a.start.cmp(&b.start).then(a.end.cmp(&b.end)));
        }
        self.sorted = true;
    }

    #[must_use]
    pub const fn is_sorted(&self) -> bool {
        self.sorted
    }

    pub fn chroms(&self) -> impl Iterator<Item = &str> {
        self.by_chrom.keys().map(String::as_str)
    }

    #[must_use]
    pub fn get(&self, chrom: &str) -> Option<&[Interval]> {
        self.by_chrom.get(chrom).map(Vec::as_slice)
    }

    pub fn iter(&self) -> impl Iterator<Item = &Interval> {
        self.by_chrom.values().flat_map(|v| v.iter())
    }

    pub fn iter_chroms(&self) -> impl Iterator<Item = (&str, &[Interval])> {
        self.by_chrom
            .iter()
            .map(|(c, v)| (c.as_str(), v.as_slice()))
    }
}

impl FromIterator<Interval> for IntervalSet {
    fn from_iter<I: IntoIterator<Item = Interval>>(items: I) -> Self {
        let mut s = Self::new();
        s.extend(items);
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(chrom: &str, start: u64, end: u64) -> Interval {
        Interval::new(chrom, start, end).unwrap()
    }

    #[test]
    fn push_and_len() {
        let mut s = IntervalSet::new();
        assert!(s.is_empty());
        s.push(iv("chr1", 100, 200));
        s.push(iv("chr2", 300, 400));
        s.push(iv("chr1", 50, 150));
        assert_eq!(s.len(), 3);
        assert_eq!(s.get("chr1").unwrap().len(), 2);
    }

    #[test]
    fn sort_orders_within_chrom() {
        let mut s = IntervalSet::from_iter([
            iv("chr1", 300, 400),
            iv("chr1", 100, 200),
            iv("chr1", 200, 300),
        ]);
        assert!(!s.is_sorted());
        s.sort();
        let v = s.get("chr1").unwrap();
        assert_eq!((v[0].start, v[1].start, v[2].start), (100, 200, 300));
        assert!(s.is_sorted());
    }

    #[test]
    fn chroms_iter_is_lexicographic() {
        let s = IntervalSet::from_iter([iv("chr10", 0, 1), iv("chr1", 0, 1), iv("chr2", 0, 1)]);
        // BTreeMap = string-lexicographic, NOT natural-numeric. "chr10"
        // sorts before "chr2". That matches bedtools' default sort order.
        let chroms: Vec<_> = s.chroms().collect();
        assert_eq!(chroms, vec!["chr1", "chr10", "chr2"]);
    }
}
