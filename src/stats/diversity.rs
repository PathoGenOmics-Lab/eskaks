//! Site-frequency-spectrum population-genetics summaries: nucleotide diversity
//! π, Watterson's θ, and Tajima's D. Pure functions of the sample size `n` and
//! the per-site derived-allele counts, so they are unit-testable in isolation.
//!
//! Conventions (Tajima 1989, *Genetics* 123:585): all three are computed over a
//! set of segregating sites in a sample of `n` sequences. π and θ are returned as
//! *totals over the region* (sum over sites), not per-site — divide by the number
//! of analysed sites to get per-site diversity. Tajima's D is dimensionless.

/// Harmonic sums a₁ = Σ_{i=1}^{n-1} 1/i and a₂ = Σ_{i=1}^{n-1} 1/i².
fn harmonic_sums(n: usize) -> (f64, f64) {
    let mut a1 = 0.0;
    let mut a2 = 0.0;
    for i in 1..n {
        let x = i as f64;
        a1 += 1.0 / x;
        a2 += 1.0 / (x * x);
    }
    (a1, a2)
}

/// Nucleotide diversity π (Tajima's estimator): the mean number of pairwise
/// differences, summed over segregating sites. Each biallelic site with derived count
/// `k` out of its own sample size `n_i` contributes `2·k·(n_i−k) / (n_i·(n_i−1))` (the
/// fraction of the C(n_i,2) pairs that differ). Per-site `n_i` lets missing genotypes
/// shrink a site's sample instead of diluting its frequency; sites with `n_i < 2`
/// contribute nothing (no pair can differ). Returned as a total over the region.
pub fn theta_pi_varn(sites: &[(usize, usize)]) -> f64 {
    let pi: f64 = sites
        .iter()
        .filter(|&&(_, ni)| ni >= 2)
        .map(|&(k, ni)| 2.0 * (k * (ni - k)) as f64 / (ni * (ni - 1)) as f64)
        .sum();
    // Rust's f64 `Sum` seeds the fold with -0.0 (to preserve a genuinely negative sum),
    // so an empty or fully-filtered set of sites returns -0.0. π is non-negative, so
    // canonicalize to +0.0; this stops the raw genome-wide summary line from printing a
    // stray "-0.000e0" and matches the "-0.0 -> 0.0" rule the file writers already apply.
    if pi == 0.0 {
        0.0
    } else {
        pi
    }
}

/// Watterson's θ_W = S / a₁, the segregating-site estimator. Returns NaN for n < 2.
pub fn theta_watterson(n: usize, s_segregating: usize) -> f64 {
    if n < 2 {
        return f64::NAN;
    }
    let (a1, _) = harmonic_sums(n);
    s_segregating as f64 / a1
}

/// Tajima's D = (π − θ_W) / sqrt(Var(π − θ_W)), the classic SFS neutrality test
/// (Tajima 1989). D < 0: excess of rare variants (purifying selection, a recent
/// sweep, or population expansion); D > 0: excess of intermediate-frequency
/// variants (balancing selection or population structure/contraction). Returns
/// NaN when it is undefined (n < 2, no segregating sites, or zero variance).
///
/// `pi` is the value from [`theta_pi`] on the same set of sites.
pub fn tajimas_d(n: usize, s_segregating: usize, pi: f64) -> f64 {
    if n < 2 || s_segregating == 0 {
        return f64::NAN;
    }
    let nf = n as f64;
    let (a1, a2) = harmonic_sums(n);
    let b1 = (nf + 1.0) / (3.0 * (nf - 1.0));
    let b2 = 2.0 * (nf * nf + nf + 3.0) / (9.0 * nf * (nf - 1.0));
    let c1 = b1 - 1.0 / a1;
    let c2 = b2 - (nf + 2.0) / (a1 * nf) + a2 / (a1 * a1);
    let e1 = c1 / a1;
    let e2 = c2 / (a1 * a1 + a2);
    let s = s_segregating as f64;
    let variance = e1 * s + e2 * s * (s - 1.0);
    if variance <= 0.0 {
        return f64::NAN;
    }
    (pi - theta_watterson(n, s_segregating)) / variance.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn harmonic_sums_are_correct() {
        // n=4: a1 = 1 + 1/2 + 1/3 = 11/6; a2 = 1 + 1/4 + 1/9 = 49/36.
        let (a1, a2) = harmonic_sums(4);
        assert!((a1 - 11.0 / 6.0).abs() < 1e-12, "a1 {}", a1);
        assert!((a2 - 49.0 / 36.0).abs() < 1e-12, "a2 {}", a2);
    }

    #[test]
    fn diversity_matches_hand_derivation_n4() {
        // n = 4 sequences, 3 segregating sites with derived counts [1, 2, 1].
        //   π   = [2·1·3 + 2·2·2 + 2·1·3] / (4·3) = (6+8+6)/12 = 20/12 = 1.666667
        //   a1  = 11/6 = 1.833333 ; θ_W = S/a1 = 3/(11/6) = 18/11 = 1.636364
        //   a2  = 49/36 = 1.361111
        //   b1  = (n+1)/(3(n−1)) = 5/9              = 0.555556
        //   b2  = 2(n²+n+3)/(9n(n−1)) = 46/108      = 0.425926
        //   c1  = b1 − 1/a1 = 5/9 − 6/11            = 0.010101
        //   c2  = b2 − (n+2)/(a1·n) + a2/a1²
        //       = 0.425926 − 6/7.333333 + 1.361111/3.361111 = 0.012703
        //   e1  = c1/a1 = 0.005510 ; e2 = c2/(a1²+a2) = 0.002690
        //   Var = e1·S + e2·S(S−1) = 0.016529 + 0.016142 = 0.032671
        //   D   = (π − θ_W)/√Var = (1.666667 − 1.636364)/0.180751 = 0.16763
        let sites = [(1usize, 4usize), (2, 4), (1, 4)];
        let pi = theta_pi_varn(&sites);
        let tw = theta_watterson(4, sites.len());
        let d = tajimas_d(4, sites.len(), pi);
        assert!((pi - 20.0 / 12.0).abs() < 1e-9, "π {}", pi);
        assert!((tw - 18.0 / 11.0).abs() < 1e-9, "θ_W {}", tw);
        assert!((d - 0.16763).abs() < 1e-4, "Tajima's D {}", d);
    }

    #[test]
    fn all_singletons_give_negative_d() {
        // Many rare (singleton) variants → π < θ_W → D < 0 (the sweep/expansion signal).
        let n = 20;
        let sites = vec![(1usize, n); 15]; // 15 singletons in a sample of n
        let pi = theta_pi_varn(&sites);
        let d = tajimas_d(n, sites.len(), pi);
        assert!(pi < theta_watterson(n, sites.len()), "π should be < θ_W for singletons");
        assert!(d < 0.0, "Tajima's D should be negative for an excess of rare variants, got {}", d);
    }

    #[test]
    fn undefined_cases_are_nan() {
        // A site with n_i < 2 contributes nothing (no pair can differ) rather than
        // making the whole sum undefined; an empty region is 0.0.
        assert_eq!(theta_pi_varn(&[(0, 1)]), 0.0);
        assert_eq!(theta_pi_varn(&[]), 0.0);
        assert!(theta_watterson(1, 3).is_nan());
        assert!(tajimas_d(10, 0, 0.0).is_nan(), "no segregating sites -> NaN");
        assert!(tajimas_d(1, 3, 1.0).is_nan(), "n<2 -> NaN");
    }

    #[test]
    fn theta_pi_varn_uses_per_site_sample_size() {
        // Uniform n_i reduces to the classic fixed-n formula: [(1,4),(2,4),(1,4)] -> 20/12.
        assert!((theta_pi_varn(&[(1, 4), (2, 4), (1, 4)]) - 20.0 / 12.0).abs() < 1e-12);
        // A site's own n_i is the denominator: 1 derived of 2 called is 2*1*1/(2*1) = 1.0.
        assert!((theta_pi_varn(&[(1, 2)]) - 1.0).abs() < 1e-12);
        // Same derived count in a larger sample dilutes the pairwise-difference fraction.
        assert!(theta_pi_varn(&[(1, 10)]) < theta_pi_varn(&[(1, 4)]));
    }

    #[test]
    fn theta_pi_varn_returns_positive_zero_not_negative_zero() {
        // Rust's f64 Sum seeds with -0.0, so an empty or fully-filtered site set would
        // otherwise yield -0.0 and print a stray "-0.000e0" in the genome-wide summary.
        for sites in [&[][..], &[(0, 1)][..], &[(0, 1), (3, 1)][..]] {
            let pi = theta_pi_varn(sites);
            assert_eq!(pi, 0.0);
            assert!(!pi.is_sign_negative(), "π must be +0.0, got a negative zero for {sites:?}");
        }
    }
}
