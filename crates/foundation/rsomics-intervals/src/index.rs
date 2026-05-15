// coitrees uses end-inclusive [first, last] intervals. We translate at the
// boundary: half-open [start, end) → coitrees::Interval::new(start, end-1, …).

use std::collections::HashMap;

use coitrees::{COITree, Interval as CoitInterval, IntervalTree};

// coitrees swaps the query-callback's node type by backend:
//   neon / avx2  → &Interval<&'a usize>     (metadata field is `&usize`)
//   basic        → &IntervalNode<usize, u32> (metadata field is `usize`)
// Tip-toe around it with a macro so both paths type-check.
#[cfg(any(target_feature = "avx2", target_feature = "neon"))]
macro_rules! meta_id {
    ($n:ident) => {
        *$n.metadata
    };
}
#[cfg(not(any(target_feature = "avx2", target_feature = "neon")))]
macro_rules! meta_id {
    ($n:ident) => {
        $n.metadata
    };
}

use crate::interval::Interval;
use crate::set::IntervalSet;

pub struct IntervalIndex {
    per_chrom: HashMap<String, COITree<usize, u32>>,
    intervals: Vec<Interval>,
}

impl IntervalIndex {
    #[must_use]
    pub fn build(set: &IntervalSet) -> Self {
        let mut intervals: Vec<Interval> = Vec::with_capacity(set.len());
        let mut by_chrom_raw: HashMap<String, Vec<CoitInterval<usize>>> = HashMap::new();
        for (chrom, ivs) in set.iter_chroms() {
            let mut nodes = Vec::with_capacity(ivs.len());
            for iv in ivs {
                let id = intervals.len();
                intervals.push(iv.clone());
                #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
                nodes.push(CoitInterval::new(iv.start as i32, (iv.end - 1) as i32, id));
            }
            by_chrom_raw.insert(chrom.to_string(), nodes);
        }
        let per_chrom = by_chrom_raw
            .into_iter()
            .map(|(chrom, nodes)| (chrom, COITree::new(&nodes)))
            .collect();
        Self {
            per_chrom,
            intervals,
        }
    }

    /// Visit every overlapping interval via callback — avoids the `Vec`
    /// allocation of [`Self::query`].
    pub fn for_each_overlap<F: FnMut(&Interval)>(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        mut f: F,
    ) {
        if start >= end {
            return;
        }
        let Some(tree) = self.per_chrom.get(chrom) else {
            return;
        };
        #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
        tree.query(start as i32, (end - 1) as i32, |node| {
            f(&self.intervals[meta_id!(node)]);
        });
    }

    #[must_use]
    pub fn query(&self, chrom: &str, start: u64, end: u64) -> Vec<&Interval> {
        if start >= end {
            return Vec::new();
        }
        let Some(tree) = self.per_chrom.get(chrom) else {
            return Vec::new();
        };
        let mut hits: Vec<&Interval> = Vec::new();
        #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
        tree.query(start as i32, (end - 1) as i32, |node| {
            hits.push(&self.intervals[meta_id!(node)]);
        });
        hits
    }

    #[must_use]
    pub fn count_overlaps(&self, chrom: &str, start: u64, end: u64) -> usize {
        if start >= end {
            return 0;
        }
        let Some(tree) = self.per_chrom.get(chrom) else {
            return 0;
        };
        #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
        tree.query_count(start as i32, (end - 1) as i32)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(chrom: &str, start: u64, end: u64) -> Interval {
        Interval::new(chrom, start, end).unwrap()
    }

    #[test]
    fn query_returns_overlapping_intervals() {
        let s = IntervalSet::from_iter([
            iv("chr1", 100, 200),
            iv("chr1", 150, 250),
            iv("chr1", 300, 400),
        ]);
        let idx = IntervalIndex::build(&s);
        let hits = idx.query("chr1", 180, 220);
        assert_eq!(hits.len(), 2);
        let pairs: Vec<_> = hits.iter().map(|h| (h.start, h.end)).collect();
        assert!(pairs.contains(&(100, 200)));
        assert!(pairs.contains(&(150, 250)));
    }

    #[test]
    fn query_other_chrom_is_empty() {
        let s = IntervalSet::from_iter([iv("chr1", 100, 200)]);
        let idx = IntervalIndex::build(&s);
        assert!(idx.query("chr2", 100, 200).is_empty());
    }

    #[test]
    fn half_open_touching_does_not_overlap() {
        // [100,200) and [200,300) share NO bases. Index must respect that.
        let s = IntervalSet::from_iter([iv("chr1", 100, 200)]);
        let idx = IntervalIndex::build(&s);
        assert_eq!(idx.count_overlaps("chr1", 200, 300), 0);
        assert_eq!(idx.count_overlaps("chr1", 199, 200), 1);
    }

    #[test]
    fn count_matches_query_len() {
        let s = IntervalSet::from_iter([
            iv("chr1", 100, 200),
            iv("chr1", 150, 250),
            iv("chr1", 300, 400),
        ]);
        let idx = IntervalIndex::build(&s);
        assert_eq!(idx.count_overlaps("chr1", 180, 220), 2);
        assert_eq!(idx.count_overlaps("chr1", 0, 50), 0);
        assert_eq!(idx.count_overlaps("chr1", 100, 400), 3);
    }

    #[test]
    fn empty_query_range_is_no_op() {
        let s = IntervalSet::from_iter([iv("chr1", 100, 200)]);
        let idx = IntervalIndex::build(&s);
        assert!(idx.query("chr1", 150, 150).is_empty());
        assert_eq!(idx.count_overlaps("chr1", 150, 150), 0);
    }

    #[test]
    fn for_each_overlap_avoids_vec() {
        let s = IntervalSet::from_iter([iv("chr1", 100, 200), iv("chr1", 150, 250)]);
        let idx = IntervalIndex::build(&s);
        let mut count = 0;
        idx.for_each_overlap("chr1", 180, 220, |_| count += 1);
        assert_eq!(count, 2);
    }
}
