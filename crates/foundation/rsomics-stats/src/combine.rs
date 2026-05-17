use statrs::distribution::{ChiSquared, ContinuousCDF, Normal};

use crate::{Result, StatsError};

// Fisher's method — -2 Σ ln(p_i) ~ χ²(2k)
pub fn fisher_combine(pvalues: &[f64]) -> Result<f64> {
    if pvalues.is_empty() {
        return Err(StatsError::Empty);
    }
    let k = pvalues.len();
    let mut stat = 0.0_f64;
    for &p in pvalues {
        if !(0.0..=1.0).contains(&p) || p.is_nan() {
            return Err(StatsError::InvalidPValue(p));
        }
        if p == 0.0 {
            return Ok(0.0);
        }
        stat += -2.0 * p.ln();
    }
    let chi2 = ChiSquared::new(2.0 * k as f64).map_err(|e| StatsError::Statrs(e.to_string()))?;
    Ok((1.0 - chi2.cdf(stat)).clamp(0.0, 1.0))
}

pub fn stouffer_combine(pvalues: &[f64], weights: Option<&[f64]>) -> Result<f64> {
    if pvalues.is_empty() {
        return Err(StatsError::Empty);
    }
    if let Some(w) = weights
        && w.len() != pvalues.len()
    {
        return Err(StatsError::SampleTooSmall {
            n: w.len(),
            required: pvalues.len(),
        });
    }
    let normal = Normal::new(0.0, 1.0).map_err(|e| StatsError::Statrs(e.to_string()))?;
    let mut numer = 0.0_f64;
    let mut denom_sq = 0.0_f64;
    for (i, &p) in pvalues.iter().enumerate() {
        if !(0.0..=1.0).contains(&p) || p.is_nan() {
            return Err(StatsError::InvalidPValue(p));
        }
        let pp = p.clamp(1e-300, 1.0 - 1e-15);
        let z = normal.inverse_cdf(1.0 - pp);
        let wi = weights.map_or(1.0, |w| w[i]);
        numer += wi * z;
        denom_sq += wi * wi;
    }
    let z_combined = numer / denom_sq.sqrt();
    Ok((1.0 - normal.cdf(z_combined)).clamp(0.0, 1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn fisher_combine_two_significant() {
        let p = fisher_combine(&[0.05, 0.05]).unwrap();
        assert!(approx(p, 0.0175, 1e-3), "p={p}");
    }

    #[test]
    fn fisher_combine_two_non_significant_stays_non_significant() {
        let p = fisher_combine(&[0.5, 0.5]).unwrap();
        assert!(p > 0.5, "p={p}");
    }

    #[test]
    fn fisher_combine_with_zero_p_returns_zero() {
        let p = fisher_combine(&[0.0, 0.5]).unwrap();
        assert_eq!(p, 0.0);
    }

    #[test]
    fn stouffer_unit_weights_two_signif() {
        let p = stouffer_combine(&[0.05, 0.05], None).unwrap();
        assert!(approx(p, 0.0100, 1e-3), "p={p}");
    }

    #[test]
    fn stouffer_weighted_matches_unweighted_when_equal() {
        let a = stouffer_combine(&[0.01, 0.04, 0.1], None).unwrap();
        let b = stouffer_combine(&[0.01, 0.04, 0.1], Some(&[1.0, 1.0, 1.0])).unwrap();
        assert!(approx(a, b, 1e-9));
    }

    #[test]
    fn empty_inputs_error() {
        assert!(matches!(fisher_combine(&[]), Err(StatsError::Empty)));
        assert!(matches!(
            stouffer_combine(&[], None),
            Err(StatsError::Empty)
        ));
    }

    #[test]
    fn invalid_p_rejected() {
        assert!(matches!(
            fisher_combine(&[0.5, 1.5]),
            Err(StatsError::InvalidPValue(_))
        ));
    }
}
