use crate::{Result, StatsError};

pub fn bonferroni_adjust(pvalues: &[f64]) -> Result<Vec<f64>> {
    validate(pvalues)?;
    let n = pvalues.len() as f64;
    Ok(pvalues.iter().map(|p| (p * n).min(1.0)).collect())
}

/// Benjamini–Hochberg FDR control. Returns adjusted p-values in the input
/// order. Walk the rank-ordered ps from highest to lowest, applying the
/// monotonic step-up correction `min(running, n*p/rank)`.
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
    for rank in (0..n).rev() {
        let i = order[rank];
        let v = (pvalues[i] * n as f64 / (rank + 1) as f64).min(1.0);
        running = running.min(v);
        adj[i] = running;
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
    }
}
