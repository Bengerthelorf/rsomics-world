#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::needless_range_loop,
    clippy::similar_names,
    clippy::many_single_char_names
)]

#[derive(Debug, Clone, Copy)]
pub struct ScoreParams {
    pub match_score: i32,
    pub mismatch: i32,
    pub gap_open: i32,   // charged on first gap byte; positive
    pub gap_extend: i32, // charged per additional byte; positive
}

impl Default for ScoreParams {
    fn default() -> Self {
        Self {
            match_score: 1,
            mismatch: -1,
            gap_open: 5,
            gap_extend: 1,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Op {
    Match,
    Mismatch,
    Insert,
    Delete,
}

#[derive(Debug, Clone)]
pub struct Alignment {
    pub score: i32,
    pub ops: Vec<Op>,
    pub a_start: usize,
    pub b_start: usize,
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum AlignError {
    #[error("empty input sequence (a={a_len}, b={b_len})")]
    Empty { a_len: usize, b_len: usize },
}

pub type Result<T> = std::result::Result<T, AlignError>;

const NEG_INF: i32 = i32::MIN / 4;

fn sub(a: u8, b: u8, p: &ScoreParams) -> i32 {
    if a.eq_ignore_ascii_case(&b) {
        p.match_score
    } else {
        p.mismatch
    }
}

/// Global alignment; Gotoh affine-gap DP (M/E/F), O(n·m).
pub fn needleman_wunsch(a: &[u8], b: &[u8], p: &ScoreParams) -> Result<Alignment> {
    if a.is_empty() || b.is_empty() {
        return Err(AlignError::Empty {
            a_len: a.len(),
            b_len: b.len(),
        });
    }
    let n = a.len();
    let m = b.len();
    let stride = m + 1;
    let mut h = vec![NEG_INF; (n + 1) * stride];
    let mut e = vec![NEG_INF; (n + 1) * stride];
    let mut f = vec![NEG_INF; (n + 1) * stride];

    h[0] = 0;
    for j in 1..=m {
        h[j] = -(p.gap_open + p.gap_extend * j as i32);
        f[j] = h[j];
    }
    for i in 1..=n {
        h[i * stride] = -(p.gap_open + p.gap_extend * i as i32);
        e[i * stride] = h[i * stride];
    }

    for i in 1..=n {
        for j in 1..=m {
            let e_open = h[(i - 1) * stride + j] - p.gap_open - p.gap_extend; // E: insert in a
            let e_ext = e[(i - 1) * stride + j] - p.gap_extend;
            e[i * stride + j] = e_open.max(e_ext);
            let f_open = h[i * stride + j - 1] - p.gap_open - p.gap_extend; // F: delete from a
            let f_ext = f[i * stride + j - 1] - p.gap_extend;
            f[i * stride + j] = f_open.max(f_ext);
            let diag = h[(i - 1) * stride + j - 1] + sub(a[i - 1], b[j - 1], p);
            h[i * stride + j] = diag.max(e[i * stride + j]).max(f[i * stride + j]);
        }
    }

    let mut ops = Vec::with_capacity(n + m);
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        if i > 0
            && j > 0
            && h[i * stride + j] == h[(i - 1) * stride + j - 1] + sub(a[i - 1], b[j - 1], p)
        {
            ops.push(if a[i - 1].eq_ignore_ascii_case(&b[j - 1]) {
                Op::Match
            } else {
                Op::Mismatch
            });
            i -= 1;
            j -= 1;
        } else if i > 0 && h[i * stride + j] == e[i * stride + j] {
            ops.push(Op::Insert);
            i -= 1;
        } else {
            ops.push(Op::Delete);
            j -= 1;
        }
    }
    ops.reverse();
    Ok(Alignment {
        score: h[n * stride + m],
        ops,
        a_start: 0,
        b_start: 0,
    })
}

/// Local alignment; Smith-Waterman affine-gap DP, score ≥ 0.
pub fn smith_waterman(a: &[u8], b: &[u8], p: &ScoreParams) -> Result<Alignment> {
    if a.is_empty() || b.is_empty() {
        return Err(AlignError::Empty {
            a_len: a.len(),
            b_len: b.len(),
        });
    }
    let n = a.len();
    let m = b.len();
    let stride = m + 1;
    let mut h = vec![0_i32; (n + 1) * stride];
    let mut e = vec![NEG_INF; (n + 1) * stride];
    let mut f = vec![NEG_INF; (n + 1) * stride];
    let mut best = 0_i32;
    let mut best_i = 0_usize;
    let mut best_j = 0_usize;

    for i in 1..=n {
        for j in 1..=m {
            let e_open = h[(i - 1) * stride + j] - p.gap_open - p.gap_extend;
            let e_ext = e[(i - 1) * stride + j] - p.gap_extend;
            e[i * stride + j] = e_open.max(e_ext);
            let f_open = h[i * stride + j - 1] - p.gap_open - p.gap_extend;
            let f_ext = f[i * stride + j - 1] - p.gap_extend;
            f[i * stride + j] = f_open.max(f_ext);
            let diag = h[(i - 1) * stride + j - 1] + sub(a[i - 1], b[j - 1], p);
            let val = diag.max(e[i * stride + j]).max(f[i * stride + j]).max(0);
            h[i * stride + j] = val;
            if val > best {
                best = val;
                best_i = i;
                best_j = j;
            }
        }
    }

    let mut ops = Vec::new();
    let (mut i, mut j) = (best_i, best_j);
    while i > 0 && j > 0 && h[i * stride + j] > 0 {
        let diag_score = h[(i - 1) * stride + j - 1] + sub(a[i - 1], b[j - 1], p);
        if h[i * stride + j] == diag_score {
            ops.push(if a[i - 1].eq_ignore_ascii_case(&b[j - 1]) {
                Op::Match
            } else {
                Op::Mismatch
            });
            i -= 1;
            j -= 1;
        } else if h[i * stride + j] == e[i * stride + j] {
            ops.push(Op::Insert);
            i -= 1;
        } else {
            ops.push(Op::Delete);
            j -= 1;
        }
    }
    ops.reverse();
    Ok(Alignment {
        score: best,
        ops,
        a_start: i,
        b_start: j,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nw_identical_seqs_score_equals_length() {
        let p = ScoreParams::default();
        let aln = needleman_wunsch(b"ACGTACGT", b"ACGTACGT", &p).unwrap();
        assert_eq!(aln.score, 8);
        assert!(aln.ops.iter().all(|o| *o == Op::Match));
    }

    #[test]
    fn nw_mismatch_subtracts_correctly() {
        let p = ScoreParams::default();
        let aln = needleman_wunsch(b"ACGT", b"ACCT", &p).unwrap();
        assert_eq!(aln.score, 2);
        let mismatches = aln.ops.iter().filter(|o| **o == Op::Mismatch).count();
        assert_eq!(mismatches, 1);
    }

    #[test]
    fn nw_gap_penalty_applied() {
        let p = ScoreParams::default();
        let aln = needleman_wunsch(b"ACGTA", b"ACGT", &p).unwrap();
        assert_eq!(aln.score, -2);
        let gaps = aln
            .ops
            .iter()
            .filter(|o| matches!(o, Op::Insert | Op::Delete))
            .count();
        assert_eq!(gaps, 1);
    }

    #[test]
    fn sw_finds_embedded_match() {
        let p = ScoreParams::default();
        let aln = smith_waterman(b"ACGT", b"TTTACGTTTT", &p).unwrap();
        assert_eq!(aln.score, 4);
        assert!(aln.ops.iter().all(|o| *o == Op::Match));
        assert_eq!(aln.a_start, 0);
        assert_eq!(aln.b_start, 3);
    }

    #[test]
    fn sw_zero_on_no_similarity() {
        let p = ScoreParams::default();
        let aln = smith_waterman(b"AAAA", b"TTTT", &p).unwrap();
        assert_eq!(aln.score, 0);
    }

    #[test]
    fn empty_input_rejected() {
        let p = ScoreParams::default();
        assert!(matches!(
            needleman_wunsch(b"", b"ACGT", &p),
            Err(AlignError::Empty { .. })
        ));
        assert!(matches!(
            smith_waterman(b"ACGT", b"", &p),
            Err(AlignError::Empty { .. })
        ));
    }

    #[test]
    fn nw_traceback_op_count_equals_alignment_length() {
        let p = ScoreParams::default();
        let aln = needleman_wunsch(b"ACGTACGT", b"ACTACGT", &p).unwrap();
        let n_ops = aln.ops.len();
        let n_aligned_a = aln
            .ops
            .iter()
            .filter(|o| matches!(o, Op::Match | Op::Mismatch | Op::Insert))
            .count();
        let n_aligned_b = aln
            .ops
            .iter()
            .filter(|o| matches!(o, Op::Match | Op::Mismatch | Op::Delete))
            .count();
        assert_eq!(n_aligned_a, 8);
        assert_eq!(n_aligned_b, 7);
        assert!(n_ops >= 8);
    }

    #[test]
    fn case_insensitive_matching() {
        let p = ScoreParams::default();
        let aln = needleman_wunsch(b"acgt", b"ACGT", &p).unwrap();
        assert_eq!(aln.score, 4);
    }
}
