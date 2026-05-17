use crate::{PopgenError, Result};

// LD r²; genotypes[hap] = [allele_A, allele_B], 0 = major / 1 = minor
pub fn r_squared(genotypes: &[[u8; 2]]) -> Result<f64> {
    if genotypes.is_empty() {
        return Err(PopgenError::Empty);
    }
    let n = genotypes.len() as f64;
    let mut a_count = 0.0;
    let mut b_count = 0.0;
    let mut ab_count = 0.0;
    for hap in genotypes {
        let a = f64::from(hap[0]);
        let b = f64::from(hap[1]);
        a_count += a;
        b_count += b;
        ab_count += a * b;
    }
    let pa = a_count / n;
    let pb = b_count / n;
    let pab = ab_count / n;
    let d = pab - pa * pb;
    let denom = pa * (1.0 - pa) * pb * (1.0 - pb);
    if denom <= 0.0 {
        // One locus is monomorphic; r² undefined. PLINK reports 0 here.
        return Ok(0.0);
    }
    Ok((d * d / denom).clamp(0.0, 1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn perfect_ld_pair() {
        let g = vec![[0, 0], [0, 0], [1, 1], [1, 1]];
        let r2 = r_squared(&g).unwrap();
        assert!(approx(r2, 1.0, 1e-9), "{r2}");
    }

    #[test]
    fn no_ld_independent_loci() {
        let g = vec![[0, 0], [0, 1], [1, 0], [1, 1]];
        let r2 = r_squared(&g).unwrap();
        assert!(approx(r2, 0.0, 1e-9), "{r2}");
    }

    #[test]
    fn monomorphic_locus_yields_zero() {
        let g = vec![[0, 0], [0, 1], [0, 1], [0, 0]];
        let r2 = r_squared(&g).unwrap();
        assert_eq!(r2, 0.0);
    }

    #[test]
    fn empty_errors() {
        assert!(matches!(r_squared(&[]), Err(PopgenError::Empty)));
    }
}
