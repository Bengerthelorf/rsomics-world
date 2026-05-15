use statrs::distribution::{ContinuousCDF, DiscreteCDF, Hypergeometric, StudentsT};

use crate::{Result, StatsError};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Alternative {
    TwoSided,
    Greater,
    Less,
}

#[derive(Debug, Clone, Copy)]
pub struct TestResult {
    pub statistic: f64,
    pub p_value: f64,
}

fn mean(xs: &[f64]) -> f64 {
    xs.iter().sum::<f64>() / xs.len() as f64
}

fn variance(xs: &[f64], mu: f64) -> f64 {
    let n = xs.len();
    if n < 2 {
        return 0.0;
    }
    xs.iter().map(|x| (x - mu).powi(2)).sum::<f64>() / (n - 1) as f64
}

pub fn welch_t(a: &[f64], b: &[f64], alt: Alternative) -> Result<TestResult> {
    if a.len() < 2 || b.len() < 2 {
        return Err(StatsError::SampleTooSmall {
            n: a.len().min(b.len()),
            required: 2,
        });
    }
    let ma = mean(a);
    let mb = mean(b);
    let va = variance(a, ma);
    let vb = variance(b, mb);
    let na = a.len() as f64;
    let nb = b.len() as f64;
    let se2 = va / na + vb / nb;
    if se2 == 0.0 {
        return Err(StatsError::ZeroVariance);
    }
    let t = (ma - mb) / se2.sqrt();
    let df = se2.powi(2) / ((va / na).powi(2) / (na - 1.0) + (vb / nb).powi(2) / (nb - 1.0)); // Welch–Satterthwaite
    let dist = StudentsT::new(0.0, 1.0, df).map_err(|e| StatsError::Statrs(e.to_string()))?;
    let p = match alt {
        Alternative::TwoSided => 2.0 * (1.0 - dist.cdf(t.abs())),
        Alternative::Greater => 1.0 - dist.cdf(t),
        Alternative::Less => dist.cdf(t),
    };
    Ok(TestResult {
        statistic: t,
        p_value: p.clamp(0.0, 1.0),
    })
}

pub fn mann_whitney_u(a: &[f64], b: &[f64], alt: Alternative) -> Result<TestResult> {
    if a.is_empty() || b.is_empty() {
        return Err(StatsError::Empty);
    }
    let n1 = a.len() as f64;
    let n2 = b.len() as f64;

    let mut pooled: Vec<(f64, u8)> = a.iter().map(|&x| (x, 0u8)).collect();
    pooled.extend(b.iter().map(|&x| (x, 1u8)));
    pooled.sort_by(|x, y| x.0.partial_cmp(&y.0).expect("NaN in input — not supported"));

    let mut ranks = vec![0.0_f64; pooled.len()];
    let mut tie_correction = 0.0_f64;
    let mut i = 0;
    while i < pooled.len() {
        let mut j = i + 1;
        while j < pooled.len() && pooled[j].0 == pooled[i].0 {
            j += 1;
        }
        let tied = (j - i) as f64;
        let mid = (i as f64 + 1.0 + j as f64) / 2.0; // mid-rank in 1-based: (i+1 + j) / 2
        for r in ranks.iter_mut().take(j).skip(i) {
            *r = mid;
        }
        if tied > 1.0 {
            tie_correction += tied.powi(3) - tied;
        }
        i = j;
    }

    let r1: f64 = ranks
        .iter()
        .zip(pooled.iter())
        .filter(|(_, p)| p.1 == 0)
        .map(|(r, _)| *r)
        .sum();
    let u1 = r1 - n1 * (n1 + 1.0) / 2.0;
    let u2 = n1 * n2 - u1;

    let mu = n1 * n2 / 2.0;
    let n = n1 + n2;
    // Normal approximation with tie correction; exact distribution for small n not implemented.
    let sigma = (n1 * n2 / 12.0 * ((n + 1.0) - tie_correction / (n * (n - 1.0)))).sqrt();
    if sigma == 0.0 {
        return Err(StatsError::ZeroVariance);
    }
    let u = u1.min(u2); // continuity-corrected below
    let z = (u - mu + 0.5) / sigma;
    let normal = statrs::distribution::Normal::new(0.0, 1.0)
        .map_err(|e| StatsError::Statrs(e.to_string()))?;
    let p = match alt {
        Alternative::TwoSided => 2.0 * normal.cdf(z.min(-z.abs())),
        Alternative::Greater => 1.0 - normal.cdf((u1 - mu - 0.5) / sigma),
        Alternative::Less => normal.cdf((u1 - mu + 0.5) / sigma),
    };
    Ok(TestResult {
        statistic: u1,
        p_value: p.clamp(0.0, 1.0),
    })
}

/// 2×2 contingency: `[[a, b], [c, d]]`. Returns the two-sided / one-sided
/// p-value via the hypergeometric tail — no approximations.
pub fn fisher_exact_2x2(a: u64, b: u64, c: u64, d: u64, alt: Alternative) -> Result<TestResult> {
    let n = a + b + c + d;
    let kk = a + c;
    let nn = a + b;
    let hyper = Hypergeometric::new(n, kk, nn).map_err(|e| StatsError::Statrs(e.to_string()))?;
    let observed = a;
    let p = match alt {
        Alternative::Less => hyper.cdf(observed),
        Alternative::Greater => 1.0 - hyper.cdf(observed.saturating_sub(1)),
        Alternative::TwoSided => {
            let observed_p = pmf(&hyper, observed);
            let mut total = 0.0;
            let lo = nn.saturating_sub(n.saturating_sub(kk));
            let hi = nn.min(kk);
            for k in lo..=hi {
                let pk = pmf(&hyper, k);
                if pk <= observed_p + 1e-12 {
                    total += pk;
                }
            }
            total
        }
    };
    let or = if b == 0 || c == 0 {
        f64::INFINITY
    } else {
        (a as f64 * d as f64) / (b as f64 * c as f64)
    };
    Ok(TestResult {
        statistic: or,
        p_value: p.clamp(0.0, 1.0),
    })
}

fn pmf(h: &Hypergeometric, k: u64) -> f64 {
    let lo = if k == 0 { 0.0 } else { h.cdf(k - 1) };
    h.cdf(k) - lo
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn welch_identical_means_gives_high_p() {
        let r = welch_t(
            &[1.0, 2.0, 3.0, 4.0],
            &[1.0, 2.0, 3.0, 4.0],
            Alternative::TwoSided,
        )
        .unwrap();
        assert!(r.p_value > 0.99, "p={}", r.p_value);
    }

    #[test]
    fn welch_known_value() {
        // df ≈ 8, t ≈ -2 → two-sided p ≈ 0.0805
        let r = welch_t(
            &[1.0, 2.0, 3.0, 4.0, 5.0],
            &[3.0, 4.0, 5.0, 6.0, 7.0],
            Alternative::TwoSided,
        )
        .unwrap();
        assert!(approx(r.p_value, 0.0805, 1e-3), "p={}", r.p_value);
    }

    #[test]
    fn welch_one_sided_is_half_two_sided_when_signed_correctly() {
        let r = welch_t(
            &[1.0, 2.0, 3.0, 4.0, 5.0],
            &[3.0, 4.0, 5.0, 6.0, 7.0],
            Alternative::Less,
        )
        .unwrap();
        assert!(approx(r.p_value, 0.0403, 1e-3), "p={}", r.p_value);
    }

    #[test]
    fn mann_whitney_separated_samples() {
        let r =
            mann_whitney_u(&[1.0, 2.0, 3.0], &[10.0, 11.0, 12.0], Alternative::TwoSided).unwrap();
        assert!(r.p_value < 0.1, "p={}", r.p_value);
    }

    #[test]
    fn mann_whitney_identical_samples_gives_p_near_1() {
        let r = mann_whitney_u(
            &[1.0, 2.0, 3.0, 4.0, 5.0],
            &[1.0, 2.0, 3.0, 4.0, 5.0],
            Alternative::TwoSided,
        )
        .unwrap();
        assert!(r.p_value > 0.9, "p={}", r.p_value);
    }

    #[test]
    fn fisher_2x2_classic_lady_tasting_tea() {
        // p = 1/70 (classic lady-tasting-tea exact result)
        let r = fisher_exact_2x2(4, 0, 0, 4, Alternative::Greater).unwrap();
        assert!(approx(r.p_value, 1.0 / 70.0, 1e-6), "p={}", r.p_value);
    }

    #[test]
    fn fisher_2x2_two_sided_independence() {
        let r = fisher_exact_2x2(5, 5, 5, 5, Alternative::TwoSided).unwrap();
        assert!(r.p_value > 0.99, "p={}", r.p_value);
    }

    #[test]
    fn welch_too_small_sample_errors() {
        assert!(matches!(
            welch_t(&[1.0], &[1.0, 2.0], Alternative::TwoSided),
            Err(StatsError::SampleTooSmall { .. })
        ));
    }
}
