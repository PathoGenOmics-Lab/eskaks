use super::*;

#[test]
fn ln_gamma_known_values() {
    assert!((ln_gamma(1.0)).abs() < 1e-9); // Γ(1)=1
    assert!((ln_gamma(2.0)).abs() < 1e-9); // Γ(2)=1
    assert!((ln_gamma(5.0) - 24.0_f64.ln()).abs() < 1e-9); // Γ(5)=4!=24
    assert!((ln_gamma(0.5) - std::f64::consts::PI.sqrt().ln()).abs() < 1e-9); // Γ(1/2)=√π
}

/// Binomial pmf, written independently of `dist.rs` (plain factorial ratio,
/// no ln_gamma) so the mid-p tests below are checked against a second implementation.
fn ref_pmf(n: u64, p0: f64) -> Vec<f64> {
    let mut row = Vec::with_capacity(n as usize + 1);
    let mut coeff = 1.0_f64; // C(n, i), updated multiplicatively
    for i in 0..=n {
        if i > 0 {
            coeff = coeff * (n - i + 1) as f64 / i as f64;
        }
        row.push(coeff * p0.powi(i as i32) * (1.0 - p0).powi((n - i) as i32));
    }
    row
}

/// The convention this code deliberately replaced: both tails count the whole atom
/// at k. Kept here only so the calibration tests can show the gap it left.
fn whole_atom_p(k: u64, n: u64, p0: f64) -> f64 {
    let row = ref_pmf(n, p0);
    let lower: f64 = row[..=k as usize].iter().sum();
    let upper: f64 = row[k as usize..].iter().sum();
    (2.0 * lower.min(upper)).min(1.0)
}

#[test]
fn binomial_two_sided_known() {
    // Mid-p: each tail counts only HALF the atom at the observed k.
    // n=10, k=8, p0=0.5: P(X>8)=11/1024, plus half of P(X=8)=45/1024 → 33.5/1024,
    // doubled = 67/1024 = 0.0654296875. The whole-atom version gave 0.109375.
    assert!((binomial_two_sided_p(8, 10, 0.5) - 0.065_429_687_5).abs() < 1e-12);
    // Symmetric centre → the two mid-p tails are both exactly 0.5, so p = 1.
    assert!((binomial_two_sided_p(5, 10, 0.5) - 1.0).abs() < 1e-12);
    // n=20, k=1, p0=0.5: 2*(P(X=0) + P(X=1)/2) = 2*(1 + 10)/2^20.
    let expected = 2.0 * 11.0 / (1u64 << 20) as f64;
    assert!((binomial_two_sided_p(1, 20, 0.5) - expected).abs() < 1e-12);
    // At an edge (k=0 or k=n) the outer tail IS the single atom, halved and doubled,
    // so mid-p collapses to the point probability itself.
    assert!((binomial_two_sided_p(0, 5, 0.3) - 0.7_f64.powi(5)).abs() < 1e-12);
    assert!((binomial_two_sided_p(4, 4, 0.75) - 0.75_f64.powi(4)).abs() < 1e-12);
    // Degenerate inputs → NaN
    assert!(binomial_two_sided_p(3, 0, 0.5).is_nan());
    assert!(binomial_two_sided_p(2, 5, 0.0).is_nan());
    assert!(binomial_two_sided_p(6, 5, 0.5).is_nan());
}

#[test]
fn binomial_mid_p_halves_exactly_one_atom() {
    // Algebraic identity that defines the convention: doubling a tail that dropped
    // half of P(X=k) is the whole-atom p-value minus one full P(X=k), whenever the
    // whole-atom value was not capped at 1.
    for &(k, n, p0) in &[(8u64, 10u64, 0.5), (1, 20, 0.5), (30, 100, 0.2), (3, 7, 0.75), (0, 4, 0.9)]
    {
        let whole = whole_atom_p(k, n, p0);
        let atom = ref_pmf(n, p0)[k as usize];
        let mid = binomial_two_sided_p(k, n, p0);
        assert!(whole < 1.0, "k={k} n={n}: this case must not be capped");
        assert!((mid - (whole - atom)).abs() < 1e-12, "k={k} n={n}: {mid} vs {}", whole - atom);
        // Strictly more powerful, and still a probability.
        assert!(mid < whole && mid > 0.0, "k={k} n={n}: mid {mid} not below whole {whole}");
    }
    // A p-value is still a probability everywhere, including at the extremes.
    for n in 1u64..=40 {
        for k in 0..=n {
            for &p0 in &[0.05_f64, 0.5, 0.73, 0.95] {
                let p = binomial_two_sided_p(k, n, p0);
                assert!((0.0..=1.0).contains(&p), "k={k} n={n} p0={p0}: p = {p}");
            }
        }
    }
    // Symmetry under p0 = 0.5 is unchanged by the halving.
    for k in 0u64..=12 {
        let a = binomial_two_sided_p(k, 12, 0.5);
        let b = binomial_two_sided_p(12 - k, 12, 0.5);
        assert!((a - b).abs() < 1e-12, "k={k}: {a} vs {b}");
    }
}

#[test]
fn binomial_mid_p_is_calibrated_where_the_whole_atom_version_was_not() {
    // Pins the calibration this convention exists for. Over every (n, p0) cell below
    // the null distribution of k is enumerated exactly (no sampling), giving the true
    // size of the test at a nominal 0.05 and the true median p-value.
    let grid_p0 = [0.5_f64, 0.6, 0.7, 0.75, 0.8];
    let (mut size_mid, mut size_whole) = (0.0_f64, 0.0_f64);
    let (mut med_mid, mut med_whole) = (0.0_f64, 0.0_f64);
    let mut cells = 0.0_f64;
    for n in 2u64..=60 {
        for &p0 in &grid_p0 {
            let row = ref_pmf(n, p0);
            let mid: Vec<f64> = (0..=n).map(|k| binomial_two_sided_p(k, n, p0)).collect();
            let whole: Vec<f64> = (0..=n).map(|k| whole_atom_p(k, n, p0)).collect();
            for (p, acc) in [(&mid, &mut size_mid), (&whole, &mut size_whole)] {
                *acc += (0..=n as usize).filter(|&k| p[k] <= 0.05).map(|k| row[k]).sum::<f64>();
            }
            // Median p under H0: walk k in increasing p until half the mass is spent.
            for (p, acc) in [(&mid, &mut med_mid), (&whole, &mut med_whole)] {
                let mut order: Vec<usize> = (0..=n as usize).collect();
                order.sort_by(|&a, &b| p[a].partial_cmp(&p[b]).unwrap());
                let mut mass = 0.0;
                for k in order {
                    mass += row[k];
                    if mass >= 0.5 {
                        *acc += p[k];
                        break;
                    }
                }
            }
            cells += 1.0;
        }
    }
    let (size_mid, size_whole) = (size_mid / cells, size_whole / cells);
    let (med_mid, med_whole) = (med_mid / cells, med_whole / cells);
    // Measured on this grid: mid-p 0.0433 / 0.5124, whole atom 0.0277 / 0.6648.
    assert!((size_mid - 0.0433).abs() < 0.002, "mid-p size drifted: {size_mid}");
    assert!((med_mid - 0.5124).abs() < 0.01, "mid-p median drifted: {med_mid}");
    // The point of the change: centred near 0.5 and spending most of the nominal
    // level, where the whole-atom convention sat at ~0.66 and ~0.028.
    assert!((med_mid - 0.5).abs() < 0.05, "mid-p median not centred: {med_mid}");
    assert!(size_mid > 0.04, "mid-p should spend most of the 0.05 level: {size_mid}");
    assert!(med_whole > 0.6, "whole-atom median should be badly off centre: {med_whole}");
    assert!(size_whole < 0.03, "whole-atom size should be far under 0.05: {size_whole}");
    // Mid-p buys that by giving up exactness: it is NOT guaranteed conservative and
    // does overshoot the nominal level on small n. Documented, not accidental.
    let overshoot = {
        let n = 6_u64;
        let p0 = 0.4;
        let row = ref_pmf(n, p0);
        (0..=n as usize).filter(|&k| binomial_two_sided_p(k as u64, n, p0) <= 0.05).map(|k| row[k]).sum::<f64>()
    };
    assert!(overshoot > 0.05, "n=6, p0=0.4 is the documented overshoot case: {overshoot}");
    assert!((overshoot - 0.0876).abs() < 0.001, "overshoot drifted: {overshoot}");
}

#[test]
fn upper_tail_is_the_whole_atom_and_matches_an_independent_pmf() {
    // Checked against `ref_pmf` (plain factorial ratio, no ln_gamma), so this is a
    // second implementation, not the same code twice.
    for &(n, p0) in &[(5u64, 0.153_323_8), (9, 0.02), (6, 0.5), (3, 0.9)] {
        let row = ref_pmf(n, p0);
        for x in 0..=n {
            let want: f64 = row[x as usize..].iter().sum();
            let got = binomial_upper_tail_p(x, n, p0);
            assert!((got - want).abs() < 1e-12, "n={n} p0={p0} x={x}: {got} vs {want}");
        }
    }

    // The whole point mass at x, NOT mid-p: the tail must exceed the mid-p version by
    // exactly half the atom. The codon scan's convention and the per-gene test's
    // convention deliberately differ, so pin the difference.
    let (n, p0, x) = (5u64, 0.153_323_8, 2usize);
    let row = ref_pmf(n, p0);
    let mid: f64 = row[x + 1..].iter().sum::<f64>() + 0.5 * row[x];
    assert!((binomial_upper_tail_p(x as u64, n, p0) - mid - 0.5 * row[x]).abs() < 1e-12);

    // The per-codon scan's own numbers: A_c = 5 possible nonsynonymous changes at a
    // TCG codon, 2 observed, theta = 140/9131 (the bundled toy genome's plug-in rate).
    let theta = 140.0 / 9131.0;
    let p = binomial_upper_tail_p(2, 5, theta);
    assert!((p - 0.002_279_558_259_877_733_6).abs() < 1e-15, "toy gene06 S44: {p}");

    // Boundaries. P(X >= 0) is 1 for every n and p0, including the empty binomial.
    assert_eq!(binomial_upper_tail_p(0, 9, 0.3), 1.0);
    assert_eq!(binomial_upper_tail_p(0, 0, 0.3), 1.0);
    assert_eq!(binomial_upper_tail_p(1, 4, 0.0), 0.0); // no change is possible
    assert_eq!(binomial_upper_tail_p(4, 4, 1.0), 1.0); // every change is certain
    assert!(binomial_upper_tail_p(5, 4, 0.3).is_nan(), "x > n is undefined");
    assert!(binomial_upper_tail_p(1, 4, f64::NAN).is_nan());
    assert!(binomial_upper_tail_p(1, 4, 1.5).is_nan());
    // Monotone in x, and a rarer null makes the same count more surprising.
    assert!(binomial_upper_tail_p(3, 9, 0.1) < binomial_upper_tail_p(2, 9, 0.1));
    assert!(binomial_upper_tail_p(2, 9, 0.01) < binomial_upper_tail_p(2, 9, 0.1));
}

#[test]
fn poisson_upper_tail_matches_a_hand_sum_and_survives_the_deep_tail() {
    // Against a direct sum of the pmf, written here without ln_gamma so the two are
    // independent implementations.
    let ref_tail = |x: u64, l: f64| -> f64 {
        let (mut term, mut acc) = ((-l).exp(), 0.0f64);
        for i in 0..2_000u64 {
            if i >= x {
                acc += term;
            }
            term = term * l / (i + 1) as f64;
        }
        acc
    };
    for l in [0.03f64, 0.5, 1.0, 3.7, 25.0] {
        for x in [1u64, 2, 3, 5, 10, 30] {
            let got = poisson_upper_tail_p(x, l);
            let want = ref_tail(x, l);
            assert!(
                (got - want).abs() <= 1e-12 + 1e-9 * want,
                "P(X>={x}) at lambda {l}: {got} vs {want}"
            );
        }
    }

    // The case the feature exists for: a codon with A_c = 5 possible nonsynonymous
    // changes, lambda = 0.0058 per change, and 5 supported independent origins of one
    // allele. `1 - CDF` would return 0 here; the upward sum keeps the value.
    let p = poisson_upper_tail_p(5, 5.0 * 0.0058);
    assert!(p > 0.0 && p < 1e-9, "five origins must be a real, tiny p: {p}");
    let deep = poisson_upper_tail_p(20, 0.03);
    assert!(deep > 0.0 && deep < 1e-45, "twenty origins at lambda 0.03: {deep}");
    assert!(deep.is_finite());

    // Boundaries and monotonicity.
    assert_eq!(poisson_upper_tail_p(0, 2.5), 1.0);
    assert_eq!(poisson_upper_tail_p(0, 0.0), 1.0);
    assert_eq!(poisson_upper_tail_p(3, 0.0), 0.0, "nothing can happen at rate 0");
    assert!(poisson_upper_tail_p(1, -1.0).is_nan());
    assert!(poisson_upper_tail_p(1, f64::NAN).is_nan());
    assert!(poisson_upper_tail_p(3, 0.1) < poisson_upper_tail_p(2, 0.1));
    assert!(poisson_upper_tail_p(2, 0.01) < poisson_upper_tail_p(2, 0.1));

    // And it agrees with the binomial it generalises: Binomial(n, p) tends to
    // Poisson(n*p) as n grows with n*p fixed, which is exactly the claim that the
    // origin null is the allele-count null one step up.
    let (n, np) = (100_000u64, 0.05f64);
    for x in [1u64, 2, 3] {
        let (a, b) = (binomial_upper_tail_p(x, n, np / n as f64), poisson_upper_tail_p(x, np));
        // The two differ only by the O(n·p²) convergence error, ~1e-5 relative here.
        assert!((a - b).abs() < 1e-4 * b.max(1e-30), "x={x}: binomial {a} vs poisson {b}");
    }
}

#[test]
fn bh_with_an_explicit_family_is_bh_when_the_family_is_complete() {
    // Same inputs, same answers: the per-gene path must be untouched by the codon
    // scan's need for a larger m.
    for p in [
        vec![0.01, 0.02, 0.03, 0.04, 0.05],
        vec![0.001, 0.5],
        vec![0.01, f64::NAN, 0.02],
    ] {
        let k = p.iter().filter(|v| v.is_finite()).count();
        let a = benjamini_hochberg(&p);
        for m in [k, 0] {
            // m == 0 is the understated-family case: it is floored at k, so an
            // under-count can never inflate significance.
            let b = benjamini_hochberg_with_m(&p, m);
            for (x, y) in a.iter().zip(&b) {
                assert!(
                    (x - y).abs() < 1e-15 || (x.is_nan() && y.is_nan()),
                    "m={m}: {x} vs {y} for {p:?}"
                );
            }
        }
    }

    // A family larger than the supplied tests scales every q by m/k. The codon scan
    // hands over only the codons carrying a SNP, while m is the whole coding genome.
    let q = benjamini_hochberg_with_m(&[0.001, 0.5], 1000);
    assert!((q[0] - 1.0).abs() < 1e-12, "0.001*1000/1 caps at 1: {}", q[0]);
    assert!((q[1] - 1.0).abs() < 1e-12, "rank 2 caps at 1: {}", q[1]);
    let q1 = benjamini_hochberg_with_m(&[1e-6, 0.5], 1000);
    assert!((q1[0] - 1e-3).abs() < 1e-15, "1e-6*1000/1: {}", q1[0]);
    // Still monotone in p, and still capped.
    let q2 = benjamini_hochberg_with_m(&[1e-9, 2e-9, 0.4], 1_317_000);
    assert!(q2[0] <= q2[1] && q2[1] <= q2[2]);
    assert!((q2[2] - 1.0).abs() < 1e-12);
}

#[test]
fn bh_uniform_and_monotone() {
    // Equal-spaced p give equal q here
    let q = benjamini_hochberg(&[0.01, 0.02, 0.03, 0.04, 0.05]);
    for v in &q {
        assert!((v - 0.05).abs() < 1e-12, "got {}", v);
    }
    // Ordering / monotonicity
    let q2 = benjamini_hochberg(&[0.001, 0.5]);
    assert!((q2[0] - 0.002).abs() < 1e-12);
    assert!((q2[1] - 0.5).abs() < 1e-12);
}

#[test]
fn bh_and_bonferroni_skip_nan() {
    let q = benjamini_hochberg(&[0.01, f64::NAN, 0.02]);
    assert!(q[1].is_nan());
    // m = 2 tested → 0.01*2/1=0.02 and 0.02*2/2=0.02
    assert!((q[0] - 0.02).abs() < 1e-12);
    assert!((q[2] - 0.02).abs() < 1e-12);

    let b = bonferroni(&[0.01, f64::NAN, 0.5]);
    assert!((b[0] - 0.02).abs() < 1e-12); // m=2
    assert!(b[1].is_nan());
    assert!((b[2] - 1.0).abs() < 1e-12); // 0.5*2 capped at 1
}

#[test]
fn rng_and_percentile() {
    // Deterministic for a fixed seed
    let mut a = SplitMix64::new(42);
    let mut b = SplitMix64::new(42);
    assert_eq!(a.next_u64(), b.next_u64());
    // below() stays in range
    let mut r = SplitMix64::new(7);
    for _ in 0..1000 {
        assert!(r.below(5) < 5);
    }
    // Percentiles
    let s = [1.0, 2.0, 3.0, 4.0, 5.0];
    assert!((percentile_sorted(&s, 0.0) - 1.0).abs() < 1e-12);
    assert!((percentile_sorted(&s, 100.0) - 5.0).abs() < 1e-12);
    assert!((percentile_sorted(&s, 50.0) - 3.0).abs() < 1e-12);
    assert!(percentile_sorted(&[], 50.0).is_nan());
}

#[test]
fn normal_p_known_values() {
    assert!((normal_two_sided_p(0.0) - 1.0).abs() < 1e-6);
    assert!((normal_two_sided_p(1.959_964) - 0.05).abs() < 1e-3);
    assert!((normal_two_sided_p(2.575_829) - 0.01).abs() < 1e-3);
    assert!(normal_two_sided_p(f64::NAN).is_nan());
    // Symmetric in sign
    assert!((normal_two_sided_p(1.5) - normal_two_sided_p(-1.5)).abs() < 1e-12);
    // Deep tail must stay tiny-but-nonzero (erfc form), not underflow to exactly 0 via
    // `1 - erf` cancellation (erf(|z|/sqrt2) rounds to 1.0 for z >~ 8.7). This is the
    // genomic-control p_gc = normal_two_sided_p(sqrt(chi2/lambda)) path.
    let p9 = normal_two_sided_p(9.0);
    assert!(p9 > 0.0 && p9 < 1e-15, "deep-tail p should be small positive, got {p9}");
    // Monotone into the tail; a more extreme Z is strictly more significant, never 0 == 0.
    assert!(normal_two_sided_p(12.0) > 0.0 && normal_two_sided_p(12.0) < p9);
}

#[test]
fn fisher_known_values() {
    // Classic tea-tasting 2x2 [[3,1],[1,3]] → two-sided p = 0.4857...
    assert!((fisher_exact_two_sided(3, 1, 1, 3) - 0.485_714_285_7).abs() < 1e-6);
    // Independent table → p near 1
    assert!((fisher_exact_two_sided(5, 5, 5, 5) - 1.0).abs() < 1e-9);
    // Strong association → small p
    assert!(fisher_exact_two_sided(10, 0, 0, 10) < 1e-3);
    // Empty table → NaN
    assert!(fisher_exact_two_sided(0, 0, 0, 0).is_nan());
    // A zero row/col is still well-defined (p ≤ 1)
    let p = fisher_exact_two_sided(4, 0, 2, 3);
    assert!((0.0..=1.0).contains(&p));
}

#[test]
fn neglog10p_matches_probability_in_representable_range() {
    // Where the raw p-value does not underflow, −log10(p) must agree. The two are the
    // same mid-p test on two scales, so this sweep is what stops one of them being
    // "fixed" without the other: every k, several n, several p0.
    for n in [1u64, 2, 5, 10, 37, 100] {
        for k in 0..=n {
            for &p0 in &[0.1_f64, 0.5, 0.62, 0.9] {
                let p = binomial_two_sided_p(k, n, p0);
                let nl = binomial_two_sided_neglog10p(k, n, p0);
                if p > 0.0 {
                    let want = -p.log10();
                    assert!(
                        (nl - want).abs() < 1e-9,
                        "k={k} n={n} p0={p0}: neglog10p {nl} vs {want} from p {p}"
                    );
                }
            }
        }
    }
    // p == 1 (the two mid-p tails balance) → −log10(p) == 0, never a tiny negative.
    assert_eq!(binomial_two_sided_neglog10p(5, 10, 0.5), 0.0);
    // Undefined cases → NaN.
    assert!(binomial_two_sided_neglog10p(0, 0, 0.5).is_nan());
    assert!(binomial_two_sided_neglog10p(11, 10, 0.5).is_nan());
}

#[test]
fn neglog10p_stays_finite_when_raw_p_underflows() {
    // 4000 nonsynonymous SNPs out of 4000 under a 0.5 null: the exact p is
    // astronomically small (underflows f64 to 0), but −log10(p) is finite.
    let p = binomial_two_sided_p(4000, 4000, 0.5);
    assert_eq!(p, 0.0, "raw p should underflow to 0");
    let nl = binomial_two_sided_neglog10p(4000, 4000, 0.5);
    assert!(nl.is_finite() && nl > 1000.0, "−log10(p) should be a large finite value, got {nl}");
    // ≈ 4000·log10(2) since P(X=4000)=0.5^4000 dominates the tail.
    assert!((nl - 4000.0 * 2.0_f64.log10()).abs() < 1.0);
}

#[test]
fn chi2_from_neglog10p_stays_finite_in_far_tail() {
    // Representable range: matches (Φ⁻¹(1−p/2))².
    for &p in &[0.5_f64, 0.05, 1e-4, 1e-8] {
        let nlp = -p.log10();
        let z = inv_normal_cdf(1.0 - p / 2.0);
        assert!((chi2_from_two_sided_neglog10p(nlp) - z * z).abs() < 1e-6);
    }
    // Far tail where 1 − p/2 rounds to 1.0: must be a large FINITE value.
    let c = chi2_from_two_sided_neglog10p(28.0); // p ≈ 1e-28
    assert!(c.is_finite() && c > 100.0 && c < 200.0, "got {c}");
    // Monotone increasing in the statistic; NaN for bad input.
    assert!(chi2_from_two_sided_neglog10p(50.0) > chi2_from_two_sided_neglog10p(20.0));
    assert!(chi2_from_two_sided_neglog10p(f64::NAN).is_nan());
    assert!(chi2_from_two_sided_neglog10p(-1.0).is_nan());
}

#[test]
fn percentile_handles_infinity_endpoints() {
    // Median lands on an ∞ element: must return ∞, not NaN (∞·0 hazard).
    assert!(percentile_sorted(&[1.0, f64::INFINITY, f64::INFINITY], 50.0).is_infinite());
    assert_eq!(percentile_sorted(&[2.0, 4.0], 0.0), 2.0);
    assert_eq!(percentile_sorted(&[2.0, 4.0], 100.0), 4.0);
    assert!((percentile_sorted(&[0.0, 10.0], 50.0) - 5.0).abs() < 1e-12);
}

#[test]
fn inv_normal_cdf_known_quantiles() {
    assert!((inv_normal_cdf(0.5)).abs() < 1e-9);
    assert!((inv_normal_cdf(0.975) - 1.959_963_984_540_054).abs() < 1e-6);
    assert!((inv_normal_cdf(0.025) + 1.959_963_984_540_054).abs() < 1e-6);
    // Symmetry and monotonicity
    assert!((inv_normal_cdf(0.8) + inv_normal_cdf(0.2)).abs() < 1e-6);
    assert!(inv_normal_cdf(0.9) > inv_normal_cdf(0.6));
    assert_eq!(inv_normal_cdf(0.0), f64::NEG_INFINITY);
    assert_eq!(inv_normal_cdf(1.0), f64::INFINITY);
}

#[test]
fn wilson_interval_properties() {
    // Contains the point estimate and is inside [0,1].
    let (lo, hi) = wilson_interval(3, 10, 0.95);
    assert!(lo >= 0.0 && hi <= 1.0 && lo < 0.3 && hi > 0.3);
    // Degenerate 0/n and n/n stay in-bounds (Wilson never leaves [0,1]).
    let (lo0, hi0) = wilson_interval(0, 8, 0.95);
    assert!(lo0 == 0.0 && hi0 > 0.0 && hi0 < 1.0);
    let (lo1, hi1) = wilson_interval(8, 8, 0.95);
    assert!(hi1 == 1.0 && lo1 > 0.0 && lo1 < 1.0);
    // n == 1 is the case that used to float one ULP below 1.0, mapping a syn == 0 gene's
    // pN/pS upper CI to a finite ~1e15 instead of +inf. Both degenerate bounds must snap.
    assert_eq!(wilson_interval(1, 1, 0.95).1, 1.0);
    assert_eq!(wilson_interval(0, 1, 0.95).0, 0.0);
    // More data → tighter interval.
    let (a, b) = wilson_interval(50, 100, 0.95);
    let (c, d) = wilson_interval(5, 10, 0.95);
    assert!((b - a) < (d - c));
    // n == 0 → NaN.
    assert!(wilson_interval(0, 0, 0.95).0.is_nan());
}

#[test]
fn cov_floataccum_default_matches_new() {
    // FloatAccum::default() delegates to new(): zeroed sums/counts, +-inf extremes.
    let d = FloatAccum::default();
    assert_eq!(d.valid_count, 0);
    assert_eq!(d.finite_ratio_count, 0);
    assert_eq!(d.sum_dn, 0.0);
    assert_eq!(d.sum_ds, 0.0);
    assert_eq!(d.sum_ratio, 0.0);
    assert_eq!(d.min_dn, f64::INFINITY);
    assert_eq!(d.max_dn, f64::NEG_INFINITY);
    assert_eq!(d.min_ratio, f64::INFINITY);
    assert_eq!(d.max_ratio, f64::NEG_INFINITY);
}

#[test]
fn cov_summarystats_default_is_empty() {
    // SummaryStats::default() delegates to new(): all counters at zero.
    let s = SummaryStats::default();
    assert_eq!(s.total_count.load(std::sync::atomic::Ordering::Relaxed), 0);
    assert_eq!(s.nan_count.load(std::sync::atomic::Ordering::Relaxed), 0);
    let f = s.floats.lock().unwrap();
    assert_eq!(f.valid_count, 0);
    assert_eq!(f.finite_ratio_count, 0);
}

#[test]
fn cov_print_summary_full_path_no_excluded() {
    // Drive print_summary through total>0, valid>0, finite ratios present, and
    // excluded == 0 (every valid pair also has a finite ratio) so the
    // `String::new()` else-branch runs, plus the full histogram loop.
    let stats = SummaryStats::new();
    stats.record_pair_atomic(0.10, 0.20, 0.5); // ratio 0.5 -> histogram bin 2
    stats.record_pair_atomic(0.30, 0.30, 1.0); // ratio 1.0 -> overflow bin 5
    let mut local = FloatAccum::new();
    local.record(0.10, 0.20, 0.5);
    local.record(0.30, 0.30, 1.0);
    stats.flush_local(&local);
    stats.print_summary(); // must not panic
    assert_eq!(stats.total_count.load(std::sync::atomic::Ordering::Relaxed), 2);
    assert_eq!(stats.nan_count.load(std::sync::atomic::Ordering::Relaxed), 0);
    let f = stats.floats.lock().unwrap();
    assert_eq!(f.valid_count, 2);
    assert_eq!(f.finite_ratio_count, 2); // => excluded == 0
}

#[test]
fn cov_percentile_single_element_returns_it() {
    // A one-element slice returns that element for any percentile.
    assert_eq!(percentile_sorted(&[42.0], 0.0), 42.0);
    assert_eq!(percentile_sorted(&[42.0], 50.0), 42.0);
    assert_eq!(percentile_sorted(&[42.0], 100.0), 42.0);
    assert_eq!(percentile_sorted(&[-3.5], 2.5), -3.5);
}

#[test]
fn cov_ln_gamma_reflection_branch_small_x() {
    // x < 0.5 uses the reflection G(x)*G(1-x) = pi / sin(pi x).
    // At x = 0.25: G(0.25)*G(0.75) = pi / sin(pi/4) = pi / (sqrt(2)/2) = pi*sqrt(2),
    //   so lnG(0.25) + lnG(0.75) = ln(pi*sqrt(2)).
    let lhs = ln_gamma(0.25) + ln_gamma(0.75);
    let rhs = (std::f64::consts::PI * std::f64::consts::SQRT_2).ln();
    assert!((lhs - rhs).abs() < 1e-9, "lhs={lhs} rhs={rhs}");
    // Known constant G(1/4) = 3.6256099082219083 => lnG(1/4) = 1.2880225246980774.
    assert!((ln_gamma(0.25) - 1.288_022_524_698_077_4).abs() < 1e-6);
}

#[test]
fn cov_inv_normal_cdf_tail_branches() {
    // Lower tail (p < 0.02425): Phi^-1(0.01) = -2.3263478740408408 (1st percentile).
    assert!((inv_normal_cdf(0.01) + 2.326_347_874_040_840_8).abs() < 1e-6);
    // Upper tail (p > 0.97575): Phi^-1(0.99) = +2.3263478740408408 (99th percentile).
    assert!((inv_normal_cdf(0.99) - 2.326_347_874_040_840_8).abs() < 1e-6);
    // Antisymmetry across the two tail branches.
    assert!((inv_normal_cdf(0.01) + inv_normal_cdf(0.99)).abs() < 1e-6);
    // Deeper into the lower tail is more negative (monotone).
    assert!(inv_normal_cdf(0.001) < inv_normal_cdf(0.01));
}

#[test]
fn cov_erfc_reflection_and_value() {
    // erfc reflects as erfc(-x) = 2 - erfc(x).
    assert!((erfc(-1.0) - (2.0 - erfc(1.0))).abs() < 1e-12);
    // Sanity vs tabulated erfc(1) = 1 - 0.8427007929 = 0.1572992071 (A&S 7.1.26).
    assert!((erfc(1.0) - 0.157_299_207_1).abs() < 1e-6);
    // Complementary to erf: 1 - erfc(1) matches the tabulated erf(1) = 0.8427007929.
    assert!(((1.0 - erfc(1.0)) - 0.842_700_792_9).abs() < 1e-6);
    // Complementary sum: erfc(x) + erfc(-x) = 2 exactly (reflection identity).
    assert!((erfc(2.5) + erfc(-2.5) - 2.0).abs() < 1e-12);
}

#[test]
fn cov_chi2_extreme_input_hits_nonfinite_break() {
    // A gigantic -log10(p) overflows ln(p) to -inf, so p underflows to 0 and the
    // far-tail fixed-point update is non-finite on the first step; the
    // `!next.is_finite()` guard breaks and the pre-loop c (= -2*ln_p = +inf)
    // is returned -- a large value, never NaN.
    let c = chi2_from_two_sided_neglog10p(1e308);
    assert!(c.is_infinite() && c > 0.0, "got {c}");
}
