use crate::{PopgenError, Result};

/// Wigginton–Cutler–Abecasis exact Hardy-Weinberg p-value for a biallelic
/// locus. `obs_hets` = observed heterozygotes, `obs_hom1`/`obs_hom2` =
/// observed homozygote counts. Returns the two-sided exact p-value.
///
/// Reference: Wigginton, Cutler & Abecasis (2005), AJHG 76:887–893.
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

    // Build the het-count distribution under HWE conditional on the observed
    // rare-allele count, then sum p(k) ≤ p(observed).
    // het parity: hets and rare allele share parity.
    let mut het_probs: Vec<f64> = vec![0.0; (rare as usize) + 1];
    // Start from the mid-point (closest valid het count).
    let mid: u64 = rare * (2 * total_genotypes - rare) / (2 * total_genotypes);
    let mid = if mid & 1 == rare & 1 { mid } else { mid + 1 };
    if mid > rare {
        return Ok(1.0);
    }
    het_probs[mid as usize] = 1.0;
    let mut sum = 1.0_f64;

    // Walk down (fewer hets).
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

    // Walk up (more hets).
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

    // Normalise so probabilities sum to 1.
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
        // 100 genotypes, AA=25, Aa=50, aa=25 — textbook HWE.
        let p = hwe_exact(50, 25, 25).unwrap();
        assert!(p > 0.99, "p={p}");
    }

    #[test]
    fn hw_extreme_homozygosity_low_p() {
        // Same allele frequencies as above (rare=100) but no heterozygotes
        // at all → strong departure from HWE.
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
        // Mild departure — return value should be a valid probability.
        let p = hwe_exact(30, 50, 20).unwrap();
        assert!((0.0..=1.0).contains(&p), "{p}");
    }
}
