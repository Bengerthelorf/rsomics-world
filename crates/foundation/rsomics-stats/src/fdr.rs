use crate::{Result, StatsError};

pub fn bonferroni_adjust(pvalues: &[f64]) -> Result<Vec<f64>> {
    validate(pvalues)?;
    let n = pvalues.len() as f64;
    Ok(pvalues.iter().map(|p| (p * n).min(1.0)).collect())
}

/// Benjamini–Hochberg FDR adjustment. Returns values in input order.
pub fn bh_adjust(pvalues: &[f64]) -> Result<Vec<f64>> {
    validate(pvalues)?;
    let n = pvalues.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| pvalues[i].partial_cmp(&pvalues[j]).expect("NaN p-value"));
    let mut adj = vec![0.0_f64; n];
    let mut running = 1.0_f64;
    for (rank, &i) in order.iter().enumerate().rev() {
        let v = (pvalues[i] * n as f64 / (rank + 1) as f64).min(1.0);
        running = running.min(v);
        adj[i] = running;
    }
    Ok(adj)
}

/// Holm step-down. R: `pmin(1, cummax((n-i+1)*p[order(p)]))[ro]`.
pub fn holm_adjust(pvalues: &[f64]) -> Result<Vec<f64>> {
    validate(pvalues)?;
    let n = pvalues.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| pvalues[i].partial_cmp(&pvalues[j]).expect("NaN p-value"));
    let mut adj = vec![0.0_f64; n];
    let mut running = 0.0_f64;
    for (rank, &i) in order.iter().enumerate() {
        let v = ((n - rank) as f64 * pvalues[i]).min(1.0);
        running = running.max(v);
        adj[i] = running;
    }
    Ok(adj)
}

/// Hochberg step-up. R: `pmin(1, cummin((n-i+1)*p[order(p,decr=T)]))[ro]`.
pub fn hochberg_adjust(pvalues: &[f64]) -> Result<Vec<f64>> {
    validate(pvalues)?;
    let n = pvalues.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| pvalues[j].partial_cmp(&pvalues[i]).expect("NaN p-value"));
    let mut adj = vec![0.0_f64; n];
    let mut running = 1.0_f64;
    for (k, &i) in order.iter().enumerate() {
        let v = ((k + 1) as f64 * pvalues[i]).min(1.0);
        running = running.min(v);
        adj[i] = running;
    }
    Ok(adj)
}

/// Benjamini–Yekutieli (BH scaled by the harmonic number, dependency-robust).
/// R: `pmin(1, cummin(q*n/i*p[order(p,decr=T)]))[ro]`, `q = sum(1/(1:n))`.
pub fn by_adjust(pvalues: &[f64]) -> Result<Vec<f64>> {
    validate(pvalues)?;
    let n = pvalues.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let q: f64 = (1..=n).map(|j| 1.0 / j as f64).sum();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| pvalues[j].partial_cmp(&pvalues[i]).expect("NaN p-value"));
    let mut adj = vec![0.0_f64; n];
    let mut running = 1.0_f64;
    for (k, &i) in order.iter().enumerate() {
        let denom = (n - k) as f64;
        let v = (q * n as f64 / denom * pvalues[i]).min(1.0);
        running = running.min(v);
        adj[i] = running;
    }
    Ok(adj)
}

/// `none`: identity (validated, already clamped to `[0,1]`). R's `"none"`.
pub fn none_adjust(pvalues: &[f64]) -> Result<Vec<f64>> {
    validate(pvalues)?;
    Ok(pvalues.to_vec())
}

/// Hommel (1988), a direct port of R `stats::p.adjust` `method = "hommel"`.
pub fn hommel_adjust(pvalues: &[f64]) -> Result<Vec<f64>> {
    validate(pvalues)?;
    let n = pvalues.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    if n == 1 {
        return Ok(vec![pvalues[0]]);
    }
    let mut o: Vec<usize> = (0..n).collect();
    o.sort_by(|&i, &j| pvalues[i].partial_cmp(&pvalues[j]).expect("NaN p-value"));
    let p: Vec<f64> = o.iter().map(|&i| pvalues[i]).collect();

    let init = (0..n)
        .map(|i| n as f64 * p[i] / (i + 1) as f64)
        .fold(f64::INFINITY, f64::min);
    let mut q = vec![init; n];
    let mut pa = vec![init; n];

    for m in (2..n).rev() {
        let i1 = n - m + 1;
        let q1 = ((n - m + 1)..n)
            .enumerate()
            .map(|(t, idx)| m as f64 * p[idx] / (t + 2) as f64)
            .fold(f64::INFINITY, f64::min);
        for idx in 0..i1 {
            q[idx] = (m as f64 * p[idx]).min(q1);
        }
        for idx in i1..n {
            q[idx] = q[i1 - 1];
        }
        for idx in 0..n {
            pa[idx] = pa[idx].max(q[idx]);
        }
    }
    let sorted_adj: Vec<f64> = (0..n).map(|i| pa[i].max(p[i])).collect();
    let mut adj = vec![0.0_f64; n];
    for (rank, &orig) in o.iter().enumerate() {
        adj[orig] = sorted_adj[rank];
    }
    Ok(adj)
}

fn validate(pvalues: &[f64]) -> Result<()> {
    for &p in pvalues {
        if !(0.0..=1.0).contains(&p) || p.is_nan() {
            return Err(StatsError::InvalidPValue(p));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_vec(a: &[f64], b: &[f64], eps: f64) -> bool {
        a.len() == b.len() && a.iter().zip(b.iter()).all(|(x, y)| (x - y).abs() < eps)
    }

    #[test]
    fn bonferroni_simple() {
        let adj = bonferroni_adjust(&[0.01, 0.04, 0.5]).unwrap();
        assert!(approx_vec(&adj, &[0.03, 0.12, 1.0], 1e-9), "{adj:?}");
    }

    #[test]
    fn bh_matches_textbook_step_up() {
        let adj = bh_adjust(&[0.01, 0.02, 0.04, 0.10, 0.50]).unwrap();
        let expected = [0.05, 0.05, 0.066_666_67, 0.125, 0.5];
        assert!(approx_vec(&adj, &expected, 1e-4), "{adj:?}");
    }

    #[test]
    fn bh_preserves_input_order() {
        let adj = bh_adjust(&[0.5, 0.01, 0.1]).unwrap();
        assert!(adj[0] > adj[2] && adj[2] > adj[1]);
    }

    #[test]
    fn bh_monotone_after_correction() {
        let ps = [0.001, 0.002, 0.003, 0.20, 0.21];
        let adj = bh_adjust(&ps).unwrap();
        let mut order: Vec<usize> = (0..adj.len()).collect();
        order.sort_by(|&i, &j| ps[i].partial_cmp(&ps[j]).unwrap());
        for w in order.windows(2) {
            assert!(
                adj[w[0]] <= adj[w[1]] + 1e-12,
                "non-monotone at {w:?}: {adj:?}"
            );
        }
    }

    #[test]
    fn invalid_p_rejected() {
        assert!(matches!(
            bh_adjust(&[0.5, 1.5]),
            Err(StatsError::InvalidPValue(_))
        ));
        assert!(matches!(
            bh_adjust(&[0.5, -0.1]),
            Err(StatsError::InvalidPValue(_))
        ));
    }

    #[test]
    fn empty_input_ok() {
        assert!(bh_adjust(&[]).unwrap().is_empty());
        assert!(bonferroni_adjust(&[]).unwrap().is_empty());
        assert!(holm_adjust(&[]).unwrap().is_empty());
        assert!(hochberg_adjust(&[]).unwrap().is_empty());
        assert!(by_adjust(&[]).unwrap().is_empty());
        assert!(hommel_adjust(&[]).unwrap().is_empty());
        assert!(none_adjust(&[]).unwrap().is_empty());
    }

    const P: [f64; 5] = [0.01, 0.02, 0.03, 0.04, 0.05];

    #[test]
    fn holm_matches_r() {
        // (n-i+1)*p sorted = [.05,.08,.09,.08,.05]; cummax = [.05,.08,.09,.09,.09]
        assert!(
            approx_vec(
                &holm_adjust(&P).unwrap(),
                &[0.05, 0.08, 0.09, 0.09, 0.09],
                1e-9
            ),
            "{:?}",
            holm_adjust(&P).unwrap()
        );
    }

    #[test]
    fn hochberg_matches_r() {
        // descending k*p_desc = [.05,.08,.09,.08,.05]; cummin → all .05
        assert!(approx_vec(
            &hochberg_adjust(&P).unwrap(),
            &[0.05, 0.05, 0.05, 0.05, 0.05],
            1e-9
        ));
    }

    #[test]
    fn by_is_bh_times_harmonic() {
        let h: f64 = (1..=5).map(|j| 1.0 / f64::from(j)).sum();
        let bh = bh_adjust(&P).unwrap();
        let by = by_adjust(&P).unwrap();
        for (b, y) in bh.iter().zip(by.iter()) {
            assert!((y - (b * h).min(1.0)).abs() < 1e-9, "bh={b} by={y} h={h}");
        }
    }

    #[test]
    fn none_is_identity() {
        assert_eq!(none_adjust(&P).unwrap(), P.to_vec());
    }

    #[test]
    fn hommel_invariants() {
        // n=1 identity; n=2 exact (pa = min(2p1, p2), then max with p).
        assert_eq!(hommel_adjust(&[0.3]).unwrap(), vec![0.3]);
        assert!(approx_vec(
            &hommel_adjust(&[0.02, 0.03]).unwrap(),
            &[0.03, 0.03],
            1e-9
        ));
        // R invariants: p ≤ hommel ≤ 1, and hommel ≤ hochberg (uniformly
        // more powerful). Order-preserving. Larger vector.
        let v = [0.001, 0.2, 0.04, 0.5, 0.005, 0.03];
        let hm = hommel_adjust(&v).unwrap();
        let hb = hochberg_adjust(&v).unwrap();
        for k in 0..v.len() {
            assert!(hm[k] >= v[k] - 1e-12 && hm[k] <= 1.0 + 1e-12, "{hm:?}");
            assert!(
                hm[k] <= hb[k] + 1e-9,
                "hommel>hochberg at {k}: {hm:?} {hb:?}"
            );
        }
    }

    #[test]
    fn family_rejects_invalid_and_preserves_order() {
        for f in [
            holm_adjust as fn(&[f64]) -> Result<Vec<f64>>,
            hochberg_adjust,
            by_adjust,
            hommel_adjust,
            none_adjust,
        ] {
            assert!(matches!(f(&[0.5, 1.5]), Err(StatsError::InvalidPValue(_))));
        }
        // order preserved: input [0.5, 0.01, 0.1] keeps adj[1] smallest.
        let a = holm_adjust(&[0.5, 0.01, 0.1]).unwrap();
        assert!(a[1] <= a[2] && a[2] <= a[0]);
    }
}
