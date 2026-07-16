use super::*;

#[test]
fn ln_gamma_known_values() {
    assert!((ln_gamma(1.0)).abs() < 1e-9); // Γ(1)=1
    assert!((ln_gamma(2.0)).abs() < 1e-9); // Γ(2)=1
    assert!((ln_gamma(5.0) - 24.0_f64.ln()).abs() < 1e-9); // Γ(5)=4!=24
    assert!((ln_gamma(0.5) - std::f64::consts::PI.sqrt().ln()).abs() < 1e-9); // Γ(1/2)=√π
}

#[test]
fn binomial_two_sided_known() {
    // n=10, k=8, p0=0.5: P(X>=8)=56/1024, doubled = 0.109375
    assert!((binomial_two_sided_p(8, 10, 0.5) - 0.109_375).abs() < 1e-9);
    // Symmetric centre → capped at 1
    assert!((binomial_two_sided_p(5, 10, 0.5) - 1.0).abs() < 1e-12);
    // n=20, k=1, p0=0.5: 2*(21/2^20)
    let expected = 2.0 * 21.0 / (1u64 << 20) as f64;
    assert!((binomial_two_sided_p(1, 20, 0.5) - expected).abs() < 1e-12);
    // Degenerate inputs → NaN
    assert!(binomial_two_sided_p(3, 0, 0.5).is_nan());
    assert!(binomial_two_sided_p(2, 5, 0.0).is_nan());
    assert!(binomial_two_sided_p(6, 5, 0.5).is_nan());
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
    // Where the raw p-value does not underflow, −log10(p) must agree.
    for &(k, n, p0) in &[(8u64, 10u64, 0.5), (2, 10, 0.5), (30, 100, 0.2), (0, 5, 0.3)] {
        let p = binomial_two_sided_p(k, n, p0);
        let nl = binomial_two_sided_neglog10p(k, n, p0);
        if p > 0.0 {
            assert!((nl - (-p.log10())).abs() < 1e-9, "k={k} n={n}: {nl} vs {}", -p.log10());
        }
    }
    // p == 1 (k exactly at the null mean) → −log10(p) == 0.
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
fn cov_erf_odd_negative_branch() {
    // erf is odd: for x < 0 the code returns -y (y = erf(|x|)); construction
    // makes erf(-a) == -erf(a) bit-for-bit for a > 0.
    assert_eq!(erf(-1.0), -erf(1.0));
    assert_eq!(erf(-2.5), -erf(2.5));
    assert!(erf(-1.0) < 0.0);
    // Sanity vs tabulated erf(1) = 0.8427007929 (A&S 7.1.26, |err| ~1.5e-7).
    assert!((erf(1.0) - 0.842_700_792_9).abs() < 1e-6);
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
