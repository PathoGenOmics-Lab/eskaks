//! Distributions and hypothesis tests: binomial, Wilson, Fisher, BH/Bonferroni, probit.

/// Lower/upper percentile of an already-ascending-sorted slice (linear interp).
pub fn percentile_sorted(sorted: &[f64], pct: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    if sorted.len() == 1 {
        return sorted[0];
    }
    let rank = (pct / 100.0) * (sorted.len() - 1) as f64;
    let lo = rank.floor() as usize;
    let hi = rank.ceil() as usize;
    let frac = rank - lo as f64;
    // Return the exact endpoint when the rank lands on one, so an ±∞ value with a
    // zero interpolation weight can't produce `∞·0 = NaN`.
    if lo == hi || frac == 0.0 {
        return sorted[lo];
    }
    if frac == 1.0 {
        return sorted[hi];
    }
    sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

// ─── Hypothesis testing helpers (dependency-free) ────────────────────────────

/// Natural log of the Gamma function via the Lanczos approximation (g = 7,
/// n = 9 coefficients). Accurate to ~1e-13 for x > 0. Used for log binomial
/// coefficients so the exact binomial test stays stable for large gene counts.
pub fn ln_gamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // Reflection: Γ(x)·Γ(1-x) = π / sin(πx)
        let pi = std::f64::consts::PI;
        (pi / (pi * x).sin()).ln() - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let t = x + G + 0.5;
        let mut a = C[0];
        for (i, &c) in C.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

/// Log of the binomial coefficient C(n, k).
fn ln_binom_coeff(n: u64, k: u64) -> f64 {
    ln_gamma(n as f64 + 1.0) - ln_gamma(k as f64 + 1.0) - ln_gamma((n - k) as f64 + 1.0)
}

/// Two-sided exact binomial-test p-value for `k` successes in `n` trials under
/// null success probability `p0`, using the "twice the smaller tail"
/// convention: p = min(1, 2·min(P(X≤k), P(X≥k))).
///
/// Returns NaN when the test is undefined (n == 0, k > n, or p0 not in (0,1)).
pub fn binomial_two_sided_p(k: u64, n: u64, p0: f64) -> f64 {
    if n == 0 || k > n || !p0.is_finite() || p0 <= 0.0 || p0 >= 1.0 {
        return f64::NAN;
    }
    let ln_p = p0.ln();
    let ln_q = (1.0 - p0).ln();
    let pmf = |i: u64| (ln_binom_coeff(n, i) + i as f64 * ln_p + (n - i) as f64 * ln_q).exp();

    let lower: f64 = (0..=k).map(pmf).sum();
    let upper: f64 = (k..=n).map(pmf).sum();
    (2.0 * lower.min(upper)).min(1.0)
}

/// Two-sided exact binomial −log10(p), computed in log space (log-sum-exp per
/// tail) so it stays finite even when the p-value is far below the ~1e-300
/// underflow floor of [`binomial_two_sided_p`]. Same "twice the smaller tail"
/// convention. Returns NaN when undefined (n == 0, k > n, p0 ∉ (0,1)); returns
/// 0.0 when p ≥ 1.
pub fn binomial_two_sided_neglog10p(k: u64, n: u64, p0: f64) -> f64 {
    if n == 0 || k > n || !p0.is_finite() || p0 <= 0.0 || p0 >= 1.0 {
        return f64::NAN;
    }
    let ln_p = p0.ln();
    let ln_q = (1.0 - p0).ln();
    let ln_pmf = |i: u64| ln_binom_coeff(n, i) + i as f64 * ln_p + (n - i) as f64 * ln_q;
    // Numerically stable sum of a tail's probabilities in log space.
    let ln_tail = |lo: u64, hi: u64| -> f64 {
        let mut mx = f64::NEG_INFINITY;
        for i in lo..=hi {
            let v = ln_pmf(i);
            if v > mx {
                mx = v;
            }
        }
        if mx == f64::NEG_INFINITY {
            return f64::NEG_INFINITY;
        }
        let s: f64 = (lo..=hi).map(|i| (ln_pmf(i) - mx).exp()).sum();
        mx + s.ln()
    };
    // Both tails include i == k, matching binomial_two_sided_p.
    let ln_two_sided = std::f64::consts::LN_2 + ln_tail(0, k).min(ln_tail(k, n));
    if ln_two_sided >= 0.0 {
        0.0
    } else {
        -ln_two_sided / std::f64::consts::LN_10
    }
}

/// Convert a two-sided −log10(p) into the matching χ²₁ statistic, staying finite
/// in the far tail. For representable p it equals (Φ⁻¹(1 − p/2))²; deeper in the
/// tail — where `1 − p/2` would round to exactly 1.0 and [`inv_normal_cdf`] would
/// return +∞ — it uses the χ²₁ upper-tail asymptotic p ≈ e^(−c/2)·√(2/(πc)),
/// solved for c by a few fixed-point steps. Returns NaN for non-finite/negative input.
pub fn chi2_from_two_sided_neglog10p(nlp: f64) -> f64 {
    if !nlp.is_finite() || nlp < 0.0 {
        return f64::NAN;
    }
    let ln_p = -nlp * std::f64::consts::LN_10; // ln(p), two-sided
    let p = ln_p.exp();
    if p > 1e-12 {
        // 1 − p/2 is safely representable here, so the exact probit is fine.
        let z = inv_normal_cdf(1.0 - p / 2.0);
        return z * z;
    }
    // ln p = −c/2 + ½·ln(2/(π c))  ⇒  c = −2 ln p + ln(2/(π c)).
    let mut c = -2.0 * ln_p;
    for _ in 0..4 {
        let next = -2.0 * ln_p + (2.0 / (std::f64::consts::PI * c)).ln();
        if !next.is_finite() {
            break;
        }
        c = next;
    }
    c
}

/// Wilson score interval for a binomial proportion `k`/`n` at confidence `conf`
/// (e.g. 0.95). Returns `(lo, hi)` clamped to [0, 1]; `(NaN, NaN)` when n == 0.
/// Robust at small n and near 0/1, unlike the normal (Wald) interval.
pub fn wilson_interval(k: u64, n: u64, conf: f64) -> (f64, f64) {
    if n == 0 {
        return (f64::NAN, f64::NAN);
    }
    // z = two-sided normal quantile for the given confidence.
    let z = inv_normal_cdf(0.5 + conf / 2.0);
    let nf = n as f64;
    let phat = k as f64 / nf;
    let z2 = z * z;
    let denom = 1.0 + z2 / nf;
    let center = (phat + z2 / (2.0 * nf)) / denom;
    let half = (z / denom) * ((phat * (1.0 - phat) / nf) + z2 / (4.0 * nf * nf)).sqrt();
    ((center - half).clamp(0.0, 1.0), (center + half).clamp(0.0, 1.0))
}

/// Inverse standard-normal CDF (probit) via Acklam's rational approximation
/// (abs error < 1.2e-9). Returns ±∞ at 0/1.
pub fn inv_normal_cdf(p: f64) -> f64 {
    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    const A: [f64; 6] = [
        -3.969_683_028_665_376e1, 2.209_460_984_245_205e2, -2.759_285_104_469_687e2,
        1.383_577_518_672_69e2, -3.066_479_806_614_716e1, 2.506_628_277_459_239e0,
    ];
    const B: [f64; 5] = [
        -5.447_609_879_822_406e1, 1.615_858_368_580_409e2, -1.556_989_798_598_866e2,
        6.680_131_188_771_972e1, -1.328_068_155_288_572e1,
    ];
    const C: [f64; 6] = [
        -7.784_894_002_430_293e-3, -3.223_964_580_411_365e-1, -2.400_758_277_161_838e0,
        -2.549_732_539_343_734e0, 4.374_664_141_464_968e0, 2.938_163_982_698_783e0,
    ];
    const D: [f64; 4] = [
        7.784_695_709_041_462e-3, 3.224_671_290_700_398e-1, 2.445_134_137_142_996e0,
        3.754_408_661_907_416e0,
    ];
    let pl = 0.024_25;
    if p < pl {
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= 1.0 - pl {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    }
}

/// Benjamini-Hochberg FDR q-values, aligned with the input. NaN inputs
/// (untested items) map to NaN and are excluded from the m used to correct.
pub fn benjamini_hochberg(pvals: &[f64]) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..pvals.len()).filter(|&i| pvals[i].is_finite()).collect();
    let m = idx.len();
    let mut q = vec![f64::NAN; pvals.len()];
    if m == 0 {
        return q;
    }
    idx.sort_by(|&a, &b| pvals[a].partial_cmp(&pvals[b]).unwrap());
    // Step-up: q_(i) = min over ranks j >= i of ( p_(j) · m / j ), capped at 1.
    let mut running_min = f64::INFINITY;
    for rank in (1..=m).rev() {
        let i = idx[rank - 1];
        running_min = running_min.min(pvals[i] * m as f64 / rank as f64);
        q[i] = running_min.min(1.0);
    }
    q
}

/// Bonferroni-corrected p-values: p·m capped at 1, where m is the number of
/// tested (finite) p-values. NaN inputs stay NaN.
pub fn bonferroni(pvals: &[f64]) -> Vec<f64> {
    let m = pvals.iter().filter(|p| p.is_finite()).count();
    pvals
        .iter()
        .map(|&p| if p.is_finite() { (p * m as f64).min(1.0) } else { f64::NAN })
        .collect()
}

/// Error function (Abramowitz & Stegun 7.1.26, max abs error ~1.5e-7).
pub(crate) fn erf(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.327_591_1 * x.abs());
    let y = 1.0
        - (((((1.061_405_429 * t - 1.453_152_027) * t) + 1.421_413_741) * t - 0.284_496_736) * t
            + 0.254_829_592)
            * t
            * (-x * x).exp();
    if x < 0.0 {
        -y
    } else {
        y
    }
}

/// Two-sided p-value for a standard-normal Z statistic: erfc(|z|/√2).
/// Used for the Nei-Gojobori analytic neutrality (Z) test.
pub fn normal_two_sided_p(z: f64) -> f64 {
    if !z.is_finite() {
        return f64::NAN;
    }
    (1.0 - erf(z.abs() / std::f64::consts::SQRT_2)).clamp(0.0, 1.0)
}

/// Two-sided Fisher's exact test p-value for the 2×2 table
/// `[[a, b], [c, d]]`, by summing the hypergeometric probabilities of all
/// tables (with the same margins) no more likely than the observed one.
/// Returns NaN for an empty table.
pub fn fisher_exact_two_sided(a: u64, b: u64, c: u64, d: u64) -> f64 {
    let r1 = a + b; // row 1 total
    let r2 = c + d; // row 2 total
    let c1 = a + c; // col 1 total
    let c2 = b + d; // col 2 total
    let n = r1 + r2;
    if n == 0 {
        return f64::NAN;
    }
    // Constant part of the log hypergeometric probability (depends on margins).
    let ln_const = ln_gamma(r1 as f64 + 1.0)
        + ln_gamma(r2 as f64 + 1.0)
        + ln_gamma(c1 as f64 + 1.0)
        + ln_gamma(c2 as f64 + 1.0)
        - ln_gamma(n as f64 + 1.0);
    // log P of the table whose top-left cell is x (margins fixed).
    let ln_p = |x: u64| -> f64 {
        let b_ = r1 - x;
        let c_ = c1 - x;
        let d_ = r2 - c_;
        ln_const
            - ln_gamma(x as f64 + 1.0)
            - ln_gamma(b_ as f64 + 1.0)
            - ln_gamma(c_ as f64 + 1.0)
            - ln_gamma(d_ as f64 + 1.0)
    };
    let p_obs = ln_p(a).exp();
    let lo = c1.saturating_sub(r2); // max(0, c1 - r2)
    let hi = c1.min(r1); // min(c1, r1)
    let mut p_sum = 0.0;
    for x in lo..=hi {
        let p = ln_p(x).exp();
        // 1e-7 tolerance so the observed table isn't excluded by rounding.
        if p <= p_obs * (1.0 + 1e-7) {
            p_sum += p;
        }
    }
    p_sum.min(1.0)
}

