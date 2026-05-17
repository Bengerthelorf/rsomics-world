use crate::{PopgenError, Result};

// exact HWE p-value — Wigginton, Cutler & Abecasis 2005, AJHG 76:887–893
pub fn hwe_exact(obs_hets: u64, obs_hom1: u64, obs_hom2: u64) -> Result<f64> {
    let n = obs_hets + obs_hom1 + obs_hom2;
    if n == 0 {
        return Err(PopgenError::Empty);
    }
    // Convention: rare allele = 2 × hom_rare + hets.
    let obs_homr = obs_hom1.min(obs_hom2);
    let rare = 2 * obs_homr + obs_hets;
    let total_genotypes = obs_hets + obs_hom1 + obs_hom2;
    let total_alleles = 2 * total_genotypes;

    if rare == 0 || rare == total_alleles {
        return Ok(1.0);
    }

    // het-count distribution under HWE given rare-allele count; two-sided p = Σ p(k) ≤ p(obs). hets and rare allele share parity
    let mut het_probs: Vec<f64> = vec![0.0; (rare as usize) + 1];
    let mid: u64 = rare * (2 * total_genotypes - rare) / (2 * total_genotypes); // HWE midpoint
    let mid = if mid & 1 == rare & 1 { mid } else { mid + 1 };
    if mid > rare {
        return Ok(1.0);
    }
    het_probs[mid as usize] = 1.0;
    let mut sum = 1.0_f64;

    let mut curr_hets = mid;
    let mut curr_homr = (rare - mid) / 2;
    let mut curr_homc = total_genotypes - curr_hets - curr_homr;
    while curr_hets >= 2 {
        let prev_prob = het_probs[curr_hets as usize];
        let new_hets = curr_hets - 2;
        let ratio = (curr_hets * (curr_hets - 1)) as f64
            / (4.0 * (curr_homr + 1) as f64 * (curr_homc + 1) as f64);
        het_probs[new_hets as usize] = prev_prob * ratio;
        sum += het_probs[new_hets as usize];
        curr_hets = new_hets;
        curr_homr += 1;
        curr_homc += 1;
    }

    let mut curr_hets = mid;
    let mut curr_homr = (rare - mid) / 2;
    let mut curr_homc = total_genotypes - curr_hets - curr_homr;
    while curr_hets <= rare.saturating_sub(2) && curr_homr > 0 && curr_homc > 0 {
        let prev_prob = het_probs[curr_hets as usize];
        let new_hets = curr_hets + 2;
        let ratio =
            (4.0 * (curr_homr as f64) * (curr_homc as f64)) / ((new_hets * (new_hets - 1)) as f64);
        het_probs[new_hets as usize] = prev_prob * ratio;
        sum += het_probs[new_hets as usize];
        curr_hets = new_hets;
        curr_homr -= 1;
        curr_homc -= 1;
    }

    let target = het_probs[obs_hets as usize] / sum;
    let p: f64 = het_probs
        .iter()
        .filter(|&&q| q / sum <= target + 1e-15)
        .map(|q| q / sum)
        .sum();
    Ok(p.clamp(0.0, 1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hw_perfect_equilibrium_high_p() {
        // AA=25, Aa=50, aa=25: textbook HWE
        let p = hwe_exact(50, 25, 25).unwrap();
        assert!(p > 0.99, "p={p}");
    }

    #[test]
    fn hw_extreme_homozygosity_low_p() {
        // same allele freq as HWE case but zero hets → strong departure
        let p = hwe_exact(0, 50, 50).unwrap();
        assert!(p < 0.001, "p={p}");
    }

    #[test]
    fn hw_monomorphic_returns_one() {
        let p = hwe_exact(0, 100, 0).unwrap();
        assert_eq!(p, 1.0);
    }

    #[test]
    fn hw_empty_errors() {
        assert!(matches!(hwe_exact(0, 0, 0), Err(PopgenError::Empty)));
    }

    #[test]
    fn hw_typical_intermediate_value() {
        let p = hwe_exact(30, 50, 20).unwrap();
        assert!((0.0..=1.0).contains(&p), "{p}");
    }
}
