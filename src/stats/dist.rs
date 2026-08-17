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
/// coefficients so the binomial test stays stable for large gene counts.
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

/// Two-sided binomial-test p-value for `k` successes in `n` trials under null
/// success probability `p0`, using the **mid-p** convention: each tail counts only
/// HALF the point mass at the observed `k`, and the smaller tail is doubled.
///
/// ```text
/// p = min(1, 2 · min( P(X<k) + ½·P(X=k) ,  P(X>k) + ½·P(X=k) ))
/// ```
///
/// # Why mid-p, and not "twice the smaller tail"
///
/// The textbook convention counts P(X=k) whole in BOTH tails. Because the binomial
/// is discrete, that version is strictly conservative: it cannot spend the nominal
/// level, and with few trials it barely spends any of it. Simulated under the exact
/// null (k ~ Binomial(n, p0), p0 drawn in [0.60, 0.85], 4e5 draws) it gives
///
/// | SNPs in the gene | median p | fraction reaching p ≤ 0.05 |
/// |------------------|----------|----------------------------|
/// | 2 to 12          | 0.83     | 1.2%                       |
/// | 10 to 60         | 0.64     | 3.1%                       |
///
/// against the 0.50 and 5% a calibrated test would give. Most genes in a bacterial
/// pN/pS scan sit in exactly that low-count range, so the whole-point-mass version
/// was discarding real signal. The same simulation under mid-p gives
///
/// | SNPs in the gene | median p | fraction reaching p ≤ 0.05 |
/// |------------------|----------|----------------------------|
/// | 2 to 12          | 0.52     | 2.6%                       |
/// | 10 to 60         | 0.51     | 4.6%                       |
///
/// centred, and roughly twice the power at the low end. Discreteness is inherent and
/// no convention lands on exactly 5% for every (n, p0); mid-p removes the systematic
/// half-atom of slack that made the old numbers so lopsided, nothing more.
///
/// # The trade: this is NOT an exact test any more
///
/// Mid-p is not guaranteed conservative and can overshoot the nominal level at small
/// n. Enumerating the exact rejection probability under H0 across n = 2..200 and
/// p0 = 0.05..0.95, at a nominal 0.05 mid-p averages 0.047 and peaks at 0.088
/// (n = 6, p0 = 0.40), whereas the whole-point-mass version averages 0.035 and never
/// exceeds 0.050. That overshoot is the accepted price of a usable discrete test
/// (Lancaster 1961; Agresti 2001) and is the right call for a genome-wide scan, but
/// it means neither the code nor the docs may describe this as an exact test.
/// Please do not "fix" it back to whole point masses: that restores the
/// miscalibration in the first table. `stats::tests` pins both tables.
///
/// Returns NaN when the test is undefined (n == 0, k > n, or p0 not in (0,1)).
pub fn binomial_two_sided_p(k: u64, n: u64, p0: f64) -> f64 {
    if n == 0 || k > n || !p0.is_finite() || p0 <= 0.0 || p0 >= 1.0 {
        return f64::NAN;
    }
    let ln_p = p0.ln();
    let ln_q = (1.0 - p0).ln();
    let pmf = |i: u64| (ln_binom_coeff(n, i) + i as f64 * ln_p + (n - i) as f64 * ln_q).exp();

    // Half the atom at k goes to each side, so lower + upper == 1 exactly and the
    // doubled minimum can never exceed 1; the cap stays as a float-rounding guard.
    let half_at_k = 0.5 * pmf(k);
    let lower: f64 = (0..k).map(pmf).sum::<f64>() + half_at_k;
    let upper: f64 = (k.saturating_add(1)..=n).map(pmf).sum::<f64>() + half_at_k;
    (2.0 * lower.min(upper)).min(1.0)
}

/// One-sided **upper-tail** binomial p-value: `P(X >= x)` for `X ~ Binomial(n, p0)`.
///
/// Unlike [`binomial_two_sided_p`] this counts the point mass at the observed `x`
/// **in full**, i.e. it is NOT the mid-p convention, and the difference is
/// deliberate. Please do not "fix" one to match the other:
///
/// * This tail is used by the per-codon recurrence scan, where `n <= 9` (a codon has
///   3 positions x 3 alternate bases) and `p0` is tiny, so the atom at `x` dominates
///   the tail. Halving it is a near-constant factor of about 2 on every codon, which
///   changes the ranking not at all and the q-values barely, whereas the error in the
///   plug-in null (a uniform per-codon mutation rate) is orders of magnitude larger.
/// * Mid-p buys calibration for a two-sided test that would otherwise spend almost
///   none of its nominal level (see [`binomial_two_sided_p`]); a one-sided tail on a
///   tiny `p0` has no such problem, and the whole atom keeps the test conservative,
///   which is the right side to err on for a genome-wide scan of about 1e6 codons.
///
/// Returns NaN when undefined (`x > n`, or `p0` outside [0, 1]). `P(X >= 0)` is 1 for
/// every `n` and `p0`, including `n == 0`.
pub fn binomial_upper_tail_p(x: u64, n: u64, p0: f64) -> f64 {
    if x > n || !p0.is_finite() || !(0.0..=1.0).contains(&p0) {
        return f64::NAN;
    }
    if x == 0 {
        return 1.0;
    }
    if p0 == 0.0 {
        return 0.0; // x >= 1 is impossible when no change can occur
    }
    if p0 == 1.0 {
        return 1.0; // X == n almost surely, and x <= n
    }
    let ln_p = p0.ln();
    let ln_q = (1.0 - p0).ln();
    // Accumulate from i == n downwards: with a small p0 the terms grow as i falls, so
    // the smallest magnitudes are added first and none is lost to rounding.
    let mut acc = 0.0f64;
    for i in (x..=n).rev() {
        acc += (ln_binom_coeff(n, i) + i as f64 * ln_p + (n - i) as f64 * ln_q).exp();
    }
    acc.clamp(0.0, 1.0)
}

/// One-sided **upper-tail** Poisson p-value: `P(X >= x)` for `X ~ Poisson(lambda)`.
///
/// The sibling of [`binomial_upper_tail_p`], for the per-codon **origin** test. There
/// the statistic is a sum of per-allele origin counts, `E_c = Σ_a h_a`, which has no
/// upper bound of 9 the way a codon's allele count does: one allele can arise fifty
/// times. Binomial(`A_c`, θ) cannot represent that, and Poisson(`A_c`·λ) is its natural
/// limit: the same null one step up, and it reduces to it when every allele arose once.
///
/// Like the binomial tail (and unlike the two-sided mid-p test) this counts the point
/// mass at the observed `x` **in full**, which keeps it conservative. Please keep the
/// two one-sided tails consistent with each other.
///
/// Summed from `x` upward rather than as `1 - CDF`, so a tail far below the f64 rounding
/// floor of the lower sum keeps its value: at λ = 0.03 and x = 20 the answer is about
/// 1e-49, which `1 - CDF` would return as exactly 0. Terms are accumulated until the
/// remainder is negligible, which is guaranteed to happen because each term is the last
/// one times `lambda / i`.
///
/// Returns NaN when undefined (`lambda` negative or non-finite). `P(X >= 0)` is 1 for
/// every λ, and a λ of exactly 0 makes any `x >= 1` impossible.
pub fn poisson_upper_tail_p(x: u64, lambda: f64) -> f64 {
    if !lambda.is_finite() || lambda < 0.0 {
        return f64::NAN;
    }
    if x == 0 {
        return 1.0;
    }
    if lambda == 0.0 {
        return 0.0;
    }
    let ln_l = lambda.ln();
    let mut acc = 0.0f64;
    let mut i = x;
    loop {
        let term = (-lambda + i as f64 * ln_l - ln_gamma(i as f64 + 1.0)).exp();
        acc += term;
        // Past the mode the terms shrink geometrically, so the remaining tail is
        // bounded by term/(1 - lambda/(i+1)) and stopping here is safe.
        if i as f64 > lambda && (term == 0.0 || term <= acc * 1e-17) {
            break;
        }
        // A pathological lambda (one far larger than any count could be) must not spin
        // forever; the tail is 1.0 to every representable digit long before here.
        if i - x > 10_000_000 {
            return 1.0;
        }
        i += 1;
    }
    acc.clamp(0.0, 1.0)
}

/// Two-sided binomial −log10(p), computed in log space (log-sum-exp per tail) so it
/// stays finite even when the p-value is far below the ~1e-300 underflow floor of
/// [`binomial_two_sided_p`]. Same **mid-p** convention as that function, and it must
/// stay that way: the two are the same test on two scales, and `stats::tests` asserts
/// they agree wherever the linear p is representable. See [`binomial_two_sided_p`]
/// for what mid-p is and why this test needs it.
///
/// Returns NaN when undefined (n == 0, k > n, p0 ∉ (0,1)); returns 0.0 when p ≥ 1.
pub fn binomial_two_sided_neglog10p(k: u64, n: u64, p0: f64) -> f64 {
    if n == 0 || k > n || !p0.is_finite() || p0 <= 0.0 || p0 >= 1.0 {
        return f64::NAN;
    }
    let ln_p = p0.ln();
    let ln_q = (1.0 - p0).ln();
    let ln_pmf = |i: u64| ln_binom_coeff(n, i) + i as f64 * ln_p + (n - i) as f64 * ln_q;
    // Mid-p weight: full mass everywhere, half (−ln 2 in log space) at the observed k.
    let ln_term = |i: u64| ln_pmf(i) - if i == k { std::f64::consts::LN_2 } else { 0.0 };
    // Numerically stable sum of a tail's probabilities in log space. Weighting the
    // atom before the log-sum-exp keeps this subtraction-free, so a tail below the
    // f64 floor loses no precision to cancellation.
    let ln_tail = |lo: u64, hi: u64| -> f64 {
        let mut mx = f64::NEG_INFINITY;
        for i in lo..=hi {
            let v = ln_term(i);
            if v > mx {
                mx = v;
            }
        }
        if mx == f64::NEG_INFINITY {
            return f64::NEG_INFINITY;
        }
        let s: f64 = (lo..=hi).map(|i| (ln_term(i) - mx).exp()).sum();
        mx + s.ln()
    };
    // Both ranges span i == k, and ln_term halves it in each, so this is
    // P(X<k) + ½P(X=k) against P(X>k) + ½P(X=k), matching binomial_two_sided_p.
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
    // k == 0 and k == n give bounds that are algebraically exactly 0.0 and 1.0; snap
    // them so float rounding of center +/- half cannot land one ULP inside (0, 1). This
    // matters for the pN/pS upper CI: a syn == 0 gene must map q -> +inf, not a finite
    // ~1e15 when hi rounds to 0.999...9 for certain n.
    let lo = if k == 0 { 0.0 } else { (center - half).clamp(0.0, 1.0) };
    let hi = if k == n { 1.0 } else { (center + half).clamp(0.0, 1.0) };
    (lo, hi)
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
    let m = pvals.iter().filter(|p| p.is_finite()).count();
    benjamini_hochberg_with_m(pvals, m)
}

/// Benjamini-Hochberg with an **explicit family size** `m`, for a scan whose family
/// is larger than the p-values it hands over.
///
/// [`benjamini_hochberg`] derives `m` from the finite entries, which is correct only
/// when every member of the family is supplied. The per-codon recurrence scan is the
/// case where it is not: it tests every codon of the coding genome but only lists the
/// codons that carry a SNP, and correcting over the listed rows alone would understate
/// the family by orders of magnitude.
///
/// **Validity condition:** the omitted `m - k` members must all have p-values no
/// smaller than every supplied one. That holds by construction for the codon scan
/// (an omitted codon has zero observed alleles, so its upper-tail p is exactly 1.0),
/// and the step-up minima for the supplied ranks are then identical to the ones a
/// full-family run would produce. `m` is floored at the number of finite p-values, so
/// an under-stated family can never inflate significance.
///
/// NaN inputs stay NaN and are outside the family, exactly as in [`benjamini_hochberg`].
pub fn benjamini_hochberg_with_m(pvals: &[f64], m: usize) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..pvals.len()).filter(|&i| pvals[i].is_finite()).collect();
    let k = idx.len();
    let mut q = vec![f64::NAN; pvals.len()];
    if k == 0 {
        return q;
    }
    // The family can never be smaller than the tests actually supplied.
    let m = m.max(k);
    idx.sort_by(|&a, &b| pvals[a].partial_cmp(&pvals[b]).unwrap());
    // Step-up: q_(i) = min over ranks j >= i of ( p_(j) · m / j ), capped at 1.
    let mut running_min = f64::INFINITY;
    for rank in (1..=k).rev() {
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

/// Complementary error function erfc(x) = 1 - erf(x), evaluated directly from the
/// Abramowitz & Stegun 7.1.26 tail form (poly(t)·exp(-x²)) so the far tail does NOT
/// collapse to 0 through `1 - erf` cancellation: erf(x) rounds to exactly 1.0 once
/// exp(-x²) drops below ~1e-16 (x ≳ 6), which would make `1 - erf` return 0.0 for a
/// still-positive tail. Same ~1.5e-7 absolute accuracy as [`erf`].
pub(crate) fn erfc(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.327_591_1 * x.abs());
    let tail = (((((1.061_405_429 * t - 1.453_152_027) * t + 1.421_413_741) * t - 0.284_496_736)
        * t
        + 0.254_829_592)
        * t)
        * (-x * x).exp();
    if x < 0.0 {
        2.0 - tail
    } else {
        tail
    }
}

/// Two-sided p-value for a standard-normal Z statistic: erfc(|z|/√2). Uses [`erfc`]
/// directly (not `1 - erf`) so a deep-tail Z (e.g. a strongly significant gene under
/// genomic control) keeps its tiny-but-nonzero p instead of underflowing to exactly 0.
pub fn normal_two_sided_p(z: f64) -> f64 {
    if !z.is_finite() {
        return f64::NAN;
    }
    erfc(z.abs() / std::f64::consts::SQRT_2).clamp(0.0, 1.0)
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

