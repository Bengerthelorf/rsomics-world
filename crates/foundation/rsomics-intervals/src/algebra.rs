use crate::index::IntervalIndex;
use crate::interval::Interval;
use crate::set::IntervalSet;

fn mk(chrom: &str, start: u64, end: u64) -> Interval {
    Interval {
        chrom: chrom.to_string(),
        start,
        end,
        strand: None,
    }
}

/// Merge overlapping or touching intervals within each chromosome.
/// Touching intervals `[a, b)` and `[b, c)` collapse into `[a, c)` — matches `bedtools merge -d 0`.
/// Strand is dropped on output.
#[must_use]
pub fn merge(input: &IntervalSet) -> IntervalSet {
    let mut work = input.clone();
    work.sort();
    let mut out = IntervalSet::new();
    for (chrom, ivs) in work.iter_chroms() {
        if ivs.is_empty() {
            continue;
        }
        let mut cur_start = ivs[0].start;
        let mut cur_end = ivs[0].end;
        for iv in &ivs[1..] {
            if iv.start <= cur_end {
                cur_end = cur_end.max(iv.end);
            } else {
                out.push(mk(chrom, cur_start, cur_end));
                cur_start = iv.start;
                cur_end = iv.end;
            }
        }
        out.push(mk(chrom, cur_start, cur_end));
    }
    out.sort();
    out
}

/// Intersect two interval sets — each overlapping pair emits an interval clipped to the overlap.
/// Mirrors `bedtools intersect` default (duplicates retained when one `a` interval hits multiple `b`).
#[must_use]
pub fn intersect(a: &IntervalSet, b: &IntervalSet) -> IntervalSet {
    // coitrees index over b: O(n_a·log n_b + output) regardless of interval length distribution.
    // An active-set sweep degrades to O(n_a·n_b) on long-a/dense-b shapes (e.g. CNV vs SNP).
    let bx = IntervalIndex::build(b);
    let mut out = IntervalSet::new();
    for (chrom, a_ivs) in a.iter_chroms() {
        for ai in a_ivs {
            bx.for_each_overlap(chrom, ai.start, ai.end, |bi| {
                let lo = ai.start.max(bi.start);
                let hi = ai.end.min(bi.end);
                if hi > lo {
                    out.push(mk(chrom, lo, hi));
                }
            });
        }
    }
    out.sort();
    out
}

/// Subtract `b`'s coverage from `a`. Mirrors `bedtools subtract` without `-A`.
#[must_use]
pub fn subtract(a: &IntervalSet, b: &IntervalSet) -> IntervalSet {
    let b_merged = merge(b);
    let mut out = IntervalSet::new();
    for (chrom, a_ivs) in a.iter_chroms() {
        let b_ivs = b_merged.get(chrom).unwrap_or(&[]);
        let mut av: Vec<&Interval> = a_ivs.iter().collect();
        av.sort_unstable_by_key(|x| (x.start, x.end));
        // b_merged is sorted + disjoint; monotone `lo` makes this O(n+m) not O(n·m).
        let mut lo = 0usize;
        for ai in &av {
            while lo < b_ivs.len() && b_ivs[lo].end <= ai.start {
                lo += 1;
            }
            for (s, e) in subtract_one(ai, &b_ivs[lo..]) {
                out.push(mk(chrom, s, e));
            }
        }
    }
    out.sort();
    out
}

fn subtract_one(a: &Interval, b_ivs: &[Interval]) -> Vec<(u64, u64)> {
    let mut cursor = a.start;
    let mut out = Vec::new();
    for b in b_ivs {
        if b.end <= cursor {
            continue;
        }
        if b.start >= a.end {
            break;
        }
        if b.start > cursor {
            out.push((cursor, b.start.min(a.end)));
        }
        cursor = cursor.max(b.end);
        if cursor >= a.end {
            break;
        }
    }
    if cursor < a.end {
        out.push((cursor, a.end));
    }
    out
}

/// Complement: uncovered regions of each chromosome, bounded by `chrom_sizes`. Mirrors `bedtools complement -g`.
#[must_use]
pub fn complement(input: &IntervalSet, chrom_sizes: &[(String, u64)]) -> IntervalSet {
    let merged = merge(input);
    let mut out = IntervalSet::new();
    for (chrom, size) in chrom_sizes {
        let ivs = merged.get(chrom).unwrap_or(&[]);
        let mut cursor: u64 = 0;
        for iv in ivs {
            if iv.start > cursor {
                out.push(mk(chrom, cursor, iv.start));
            }
            cursor = cursor.max(iv.end);
            if cursor >= *size {
                break;
            }
        }
        if cursor < *size {
            out.push(mk(chrom, cursor, *size));
        }
    }
    out.sort();
    out
}

#[must_use]
pub fn coverage_bases(input: &IntervalSet) -> Vec<(String, u64)> {
    let merged = merge(input);
    let mut out = Vec::new();
    for (chrom, ivs) in merged.iter_chroms() {
        let total: u64 = ivs.iter().map(Interval::len).sum();
        out.push((chrom.to_string(), total));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(chrom: &str, start: u64, end: u64) -> Interval {
        Interval::new(chrom, start, end).unwrap()
    }

    fn set(items: impl IntoIterator<Item = Interval>) -> IntervalSet {
        let mut s: IntervalSet = items.into_iter().collect();
        s.sort();
        s
    }

    #[test]
    fn merge_collapses_overlapping_and_touching() {
        let s = set([
            iv("chr1", 100, 200),
            iv("chr1", 150, 250),
            iv("chr1", 250, 300),
            iv("chr1", 400, 500),
        ]);
        let m = merge(&s);
        let v = m.get("chr1").unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!((v[0].start, v[0].end), (100, 300));
        assert_eq!((v[1].start, v[1].end), (400, 500));
    }

    #[test]
    fn merge_independent_per_chrom() {
        let s = set([iv("chr1", 100, 200), iv("chr2", 100, 200)]);
        let m = merge(&s);
        assert_eq!(m.get("chr1").unwrap().len(), 1);
        assert_eq!(m.get("chr2").unwrap().len(), 1);
    }

    #[test]
    fn intersect_clips_to_overlap() {
        let a = set([iv("chr1", 100, 200), iv("chr1", 300, 400)]);
        let b = set([iv("chr1", 150, 350)]);
        let i = intersect(&a, &b);
        let v = i.get("chr1").unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!((v[0].start, v[0].end), (150, 200));
        assert_eq!((v[1].start, v[1].end), (300, 350));
    }

    #[test]
    fn intersect_long_a_then_short_a_dense_b() {
        let a = set([iv("chr1", 0, 100), iv("chr1", 10, 11)]);
        let b = set([iv("chr1", 5, 15), iv("chr1", 50, 60)]);
        let i = intersect(&a, &b);
        let v = i.get("chr1").unwrap();
        assert_eq!(
            v.iter().map(|x| (x.start, x.end)).collect::<Vec<_>>(),
            vec![(5, 15), (10, 11), (50, 60)]
        );
    }

    #[test]
    fn intersect_different_chroms_yield_nothing() {
        let a = set([iv("chr1", 100, 200)]);
        let b = set([iv("chr2", 100, 200)]);
        assert!(intersect(&a, &b).is_empty());
    }

    #[test]
    fn subtract_punches_holes() {
        let a = set([iv("chr1", 100, 500)]);
        let b = set([iv("chr1", 200, 300), iv("chr1", 400, 450)]);
        let s = subtract(&a, &b);
        let v = s.get("chr1").unwrap();
        assert_eq!(v.len(), 3);
        assert_eq!((v[0].start, v[0].end), (100, 200));
        assert_eq!((v[1].start, v[1].end), (300, 400));
        assert_eq!((v[2].start, v[2].end), (450, 500));
    }

    #[test]
    fn subtract_full_cover_yields_nothing() {
        let a = set([iv("chr1", 100, 200)]);
        let b = set([iv("chr1", 50, 250)]);
        assert!(subtract(&a, &b).is_empty());
    }

    #[test]
    fn complement_against_genome_sizes() {
        let s = set([iv("chr1", 100, 200), iv("chr1", 400, 500)]);
        let sizes = vec![("chr1".to_string(), 1000)];
        let c = complement(&s, &sizes);
        let v = c.get("chr1").unwrap();
        assert_eq!(v.len(), 3);
        assert_eq!((v[0].start, v[0].end), (0, 100));
        assert_eq!((v[1].start, v[1].end), (200, 400));
        assert_eq!((v[2].start, v[2].end), (500, 1000));
    }

    #[test]
    fn complement_chrom_with_no_intervals_is_full_length() {
        let s = IntervalSet::new();
        let sizes = vec![("chr1".to_string(), 1000)];
        let c = complement(&s, &sizes);
        let v = c.get("chr1").unwrap();
        assert_eq!(v.len(), 1);
        assert_eq!((v[0].start, v[0].end), (0, 1000));
    }

    #[test]
    fn coverage_bases_dedupes_overlap() {
        let s = set([iv("chr1", 100, 200), iv("chr1", 150, 250)]);
        let cov = coverage_bases(&s);
        assert_eq!(cov, vec![("chr1".to_string(), 150)]);
    }
}
