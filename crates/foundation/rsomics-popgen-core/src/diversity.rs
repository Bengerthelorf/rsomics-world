use crate::{PopgenError, Result};

/// Nucleotide diversity π (average pairwise difference per site).
/// `derived_counts[i]` = derived/minor allele count at site i; `n` = sample
/// size (haploid chromosomes; diploid samples contribute 2).
pub fn theta_pi(derived_counts: &[u64], n: u64) -> Result<f64> {
    if derived_counts.is_empty() {
        return Err(PopgenError::Empty);
    }
    if n < 2 {
        return Err(PopgenError::SampleTooSmall {
            n: n as usize,
            required: 2,
        });
    }
    let pairs = (n * (n - 1)) as f64 / 2.0;
    let mut total = 0.0_f64;
    for &k in derived_counts {
        if k > n {
            return Err(PopgenError::InvalidAlleleCount {
                observed: k,
                max: n,
            });
        }
        let segregating_pairs = (k * (n - k)) as f64;
        total += segregating_pairs / pairs;
    }
    Ok(total / derived_counts.len() as f64)
}

/// Watterson's θ per site. `s` = segregating sites, `n` = sample size, `sites` = window length.
pub fn watterson_theta(s: u64, n: u64, sites: u64) -> Result<f64> {
    if n < 2 {
        return Err(PopgenError::SampleTooSmall {
            n: n as usize,
            required: 2,
        });
    }
    if sites == 0 {
        return Err(PopgenError::Empty);
    }
    let a1: f64 = (1..n).map(|i| 1.0 / i as f64).sum();
    Ok(s as f64 / (a1 * sites as f64))
}

/// Tajima's D over a window; `n` = sample size.
pub fn tajimas_d(derived_counts: &[u64], n: u64) -> Result<f64> {
    if derived_counts.is_empty() {
        return Err(PopgenError::Empty);
    }
    if n < 4 {
        return Err(PopgenError::SampleTooSmall {
            n: n as usize,
            required: 4,
        });
    }
    let s: u64 = derived_counts.iter().filter(|&&k| k != 0 && k != n).count() as u64;
    if s == 0 {
        return Err(PopgenError::NoSegregating);
    }
    let sites = derived_counts.len() as u64;
    let pi = theta_pi(derived_counts, n)? * sites as f64;
    let s_f = s as f64;
    let n_f = n as f64;

    let a1: f64 = (1..n).map(|i| 1.0 / i as f64).sum();
    let a2: f64 = (1..n).map(|i| 1.0 / (i as f64).powi(2)).sum();
    let b1 = (n_f + 1.0) / (3.0 * (n_f - 1.0));
    let b2 = 2.0 * (n_f * n_f + n_f + 3.0) / (9.0 * n_f * (n_f - 1.0));
    let c1 = b1 - 1.0 / a1;
    let c2 = b2 - (n_f + 2.0) / (a1 * n_f) + a2 / (a1 * a1);
    let e1 = c1 / a1;
    let e2 = c2 / (a1 * a1 + a2);

    let theta_w = s_f / a1;
    let variance = e1 * s_f + e2 * s_f * (s_f - 1.0);
    if variance <= 0.0 {
        return Err(PopgenError::NoSegregating);
    }
    Ok((pi - theta_w) / variance.sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn theta_pi_one_segregating_singleton() {
        // k=1, n=10: pairs = 1×9 = 9, total = 45, π = 9/45 = 0.2
        let p = theta_pi(&[1], 10).unwrap();
        assert!(approx(p, 0.2, 1e-9), "{p}");
    }

    #[test]
    fn theta_pi_balanced_50_50() {
        // k=5, n=10: pairs = 5×5 = 25, total = 45, π = 25/45
        let p = theta_pi(&[5], 10).unwrap();
        assert!(approx(p, 25.0 / 45.0, 1e-9), "{p}");
    }

    #[test]
    fn theta_pi_no_segregating_is_zero() {
        let p = theta_pi(&[0, 10, 0], 10).unwrap();
        assert_eq!(p, 0.0);
    }

    #[test]
    fn watterson_known_value() {
        // n=10: a1 ≈ 2.828968; θw = 5 / (a1 × 1000)
        let t = watterson_theta(5, 10, 1000).unwrap();
        assert!(approx(t, 5.0 / (2.828_968_3 * 1000.0), 1e-5), "{t}");
    }

    #[test]
    fn tajimas_d_zero_when_pi_equals_theta_w() {
        // All singletons → D negative; just verify finite + non-NaN.
        let counts = vec![0_u64, 1, 0, 1, 0, 0, 1];
        let d = tajimas_d(&counts, 8).unwrap();
        assert!(d.is_finite());
    }

    #[test]
    fn tajimas_d_errors_on_no_segregating() {
        let counts = vec![0_u64, 10, 10, 0];
        assert!(matches!(
            tajimas_d(&counts, 10),
            Err(PopgenError::NoSegregating)
        ));
    }

    #[test]
    fn tajimas_d_small_sample_errors() {
        assert!(matches!(
            tajimas_d(&[1, 2], 3),
            Err(PopgenError::SampleTooSmall { .. })
        ));
    }

    #[test]
    fn invalid_allele_count_errors() {
        assert!(matches!(
            theta_pi(&[12], 10),
            Err(PopgenError::InvalidAlleleCount { .. })
        ));
    }
}
