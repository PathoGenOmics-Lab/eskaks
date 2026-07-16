use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

/// Histogram bin edges for dN/dS ratio distribution.
const BIN_EDGES: [f64; 5] = [0.2, 0.4, 0.6, 0.8, 1.0];
const BIN_LABELS: [&str; 6] = [
    "[0.0, 0.2)",
    "[0.2, 0.4)",
    "[0.4, 0.6)",
    "[0.6, 0.8)",
    "[0.8, 1.0)",
    "[1.0, inf)",
];
const NUM_BINS: usize = 6;
const MAX_BAR_WIDTH: usize = 20;

/// Thread-safe accumulator for dN/dS summary statistics.
/// Uses AtomicUsize for counters and Mutex<FloatAccum> for float aggregates.
pub struct SummaryStats {
    pub total_count: AtomicUsize,
    pub nan_count: AtomicUsize,
    pub histogram: [AtomicUsize; NUM_BINS],
    pub floats: Mutex<FloatAccum>,
}

/// Float accumulator for min/max/sum of dN, dS, and ratio.
/// Merged periodically from thread-local copies.
#[derive(Clone)]
pub struct FloatAccum {
    pub sum_dn: f64,
    pub sum_ds: f64,
    pub sum_ratio: f64,
    pub min_dn: f64,
    pub max_dn: f64,
    pub min_ds: f64,
    pub max_ds: f64,
    pub min_ratio: f64,
    pub max_ratio: f64,
    pub finite_ratio_count: usize,
    pub valid_count: usize,
}

impl Default for FloatAccum {
    fn default() -> Self { Self::new() }
}

impl FloatAccum {
    pub fn new() -> Self {
        FloatAccum {
            sum_dn: 0.0,
            sum_ds: 0.0,
            sum_ratio: 0.0,
            min_dn: f64::INFINITY,
            max_dn: f64::NEG_INFINITY,
            min_ds: f64::INFINITY,
            max_ds: f64::NEG_INFINITY,
            min_ratio: f64::INFINITY,
            max_ratio: f64::NEG_INFINITY,
            finite_ratio_count: 0,
            valid_count: 0,
        }
    }

    /// Record a single pair result into this local accumulator.
    #[inline]
    pub fn record(&mut self, dn: f64, ds: f64, ratio: f64) {
        if dn.is_finite() && ds.is_finite() {
            self.sum_dn += dn;
            self.sum_ds += ds;
            if dn < self.min_dn { self.min_dn = dn; }
            if dn > self.max_dn { self.max_dn = dn; }
            if ds < self.min_ds { self.min_ds = ds; }
            if ds > self.max_ds { self.max_ds = ds; }
            self.valid_count += 1;
        }
        if ratio.is_finite() {
            self.sum_ratio += ratio;
            if ratio < self.min_ratio { self.min_ratio = ratio; }
            if ratio > self.max_ratio { self.max_ratio = ratio; }
            self.finite_ratio_count += 1;
        }
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: &FloatAccum) {
        self.sum_dn += other.sum_dn;
        self.sum_ds += other.sum_ds;
        self.sum_ratio += other.sum_ratio;
        if other.min_dn < self.min_dn { self.min_dn = other.min_dn; }
        if other.max_dn > self.max_dn { self.max_dn = other.max_dn; }
        if other.min_ds < self.min_ds { self.min_ds = other.min_ds; }
        if other.max_ds > self.max_ds { self.max_ds = other.max_ds; }
        if other.min_ratio < self.min_ratio { self.min_ratio = other.min_ratio; }
        if other.max_ratio > self.max_ratio { self.max_ratio = other.max_ratio; }
        self.finite_ratio_count += other.finite_ratio_count;
        self.valid_count += other.valid_count;
    }

    /// Reset to initial state for reuse in next row.
    pub fn reset(&mut self) {
        *self = FloatAccum::new();
    }
}

impl Default for SummaryStats {
    fn default() -> Self { Self::new() }
}

impl SummaryStats {
    pub fn new() -> Self {
        SummaryStats {
            total_count: AtomicUsize::new(0),
            nan_count: AtomicUsize::new(0),
            histogram: [
                AtomicUsize::new(0), AtomicUsize::new(0), AtomicUsize::new(0),
                AtomicUsize::new(0), AtomicUsize::new(0), AtomicUsize::new(0),
            ],
            floats: Mutex::new(FloatAccum::new()),
        }
    }

    /// Determine histogram bin index for a dN/dS ratio.
    #[inline]
    pub fn bin_index(ratio: f64) -> usize {
        if !ratio.is_finite() || ratio < 0.0 {
            return NUM_BINS - 1; // overflow bin
        }
        for (i, &edge) in BIN_EDGES.iter().enumerate() {
            if ratio < edge {
                return i;
            }
        }
        NUM_BINS - 1
    }

    /// Record a single pair into the atomic counters (call from any thread).
    #[inline]
    pub fn record_pair_atomic(&self, dn: f64, ds: f64, ratio: f64) {
        self.total_count.fetch_add(1, Ordering::Relaxed);
        if !dn.is_finite() || !ds.is_finite() {
            self.nan_count.fetch_add(1, Ordering::Relaxed);
        } else {
            let bin = Self::bin_index(ratio);
            self.histogram[bin].fetch_add(1, Ordering::Relaxed);
        }
    }

    /// Flush a thread-local FloatAccum into the shared Mutex.
    pub fn flush_local(&self, local: &FloatAccum) {
        let mut shared = self.floats.lock().expect("SummaryStats mutex poisoned");
        shared.merge(local);
    }

    /// Print summary to stderr.
    pub fn print_summary(&self) {
        let total = self.total_count.load(Ordering::Relaxed);
        let nans = self.nan_count.load(Ordering::Relaxed);
        let valid = total.saturating_sub(nans);
        let floats = self.floats.lock().expect("SummaryStats mutex poisoned");

        eprintln!();
        eprintln!("\u{2550}\u{2550}\u{2550} dN/dS Summary \u{2550}\u{2550}\u{2550}");
        eprintln!("  Total pairs:     {}", total);
        if total > 0 {
            eprintln!("  Valid pairs:     {}  ({:.1}%)", valid, 100.0 * valid as f64 / total as f64);
            eprintln!("  Saturated (NaN): {}  ({:.1}%)", nans, 100.0 * nans as f64 / total as f64);
        }

        if floats.valid_count > 0 {
            let mean_dn = floats.sum_dn / floats.valid_count as f64;
            let mean_ds = floats.sum_ds / floats.valid_count as f64;
            eprintln!();
            eprintln!("  dN:    min={:.6}  max={:.6}  mean={:.6}", floats.min_dn, floats.max_dn, mean_dn);
            eprintln!("  dS:    min={:.6}  max={:.6}  mean={:.6}", floats.min_ds, floats.max_ds, mean_ds);
        }
        // Pooled ratio = mean(dN)/mean(dS). This is the preferred summary: the
        // mean of per-pair ratios is statistically biased (mean of ratios !=
        // ratio of means) and is dominated by pairs with a near-zero dS.
        if floats.valid_count > 0 && floats.sum_ds > 0.0 {
            eprintln!("  dN/dS (pooled = meanDN/meanDS): {:.6}", floats.sum_dn / floats.sum_ds);
        }
        if floats.finite_ratio_count > 0 {
            let mean_ratio = floats.sum_ratio / floats.finite_ratio_count as f64;
            let excluded = floats.valid_count.saturating_sub(floats.finite_ratio_count);
            let note = if excluded > 0 {
                format!(", {} undefined excluded", excluded)
            } else {
                String::new()
            };
            eprintln!(
                "  dN/dS (mean of per-pair ratios): min={:.6}  max={:.6}  mean={:.6}  (n={}{})",
                floats.min_ratio, floats.max_ratio, mean_ratio, floats.finite_ratio_count, note
            );
        }

        // Histogram
        let bins: Vec<usize> = self.histogram.iter().map(|a| a.load(Ordering::Relaxed)).collect();
        let max_count = *bins.iter().max().unwrap_or(&1);

        if valid > 0 {
            eprintln!();
            eprintln!("  dN/dS Distribution:");
            for (i, &count) in bins.iter().enumerate() {
                let bar_len = if max_count > 0 {
                    (count * MAX_BAR_WIDTH) / max_count
                } else {
                    0
                };
                let bar: String = "\u{2588}".repeat(bar_len);
                let pad: String = " ".repeat(MAX_BAR_WIDTH - bar_len);
                let pct = 100.0 * count as f64 / valid as f64;
                eprintln!("  {} {}{} {:>3} ({:>5.1}%)", BIN_LABELS[i], bar, pad, count, pct);
            }
        }

        eprintln!("\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}\u{2550}");
        eprintln!();
    }

    /// Get histogram bin counts (for SVG plotting).
    pub fn get_histogram(&self) -> Vec<(String, usize)> {
        BIN_LABELS.iter().zip(self.histogram.iter())
            .map(|(&label, count)| (label.to_string(), count.load(Ordering::Relaxed)))
            .collect()
    }
}

/// Per-window accumulator for sliding window plot data.
pub struct WindowStats {
    pub bins: Vec<Mutex<WindowBin>>,
}

#[derive(Clone)]
pub struct WindowBin {
    pub sum_ratio: f64,
    pub sum_sq: f64,
    pub count: usize,
}

impl WindowStats {
    pub fn new(num_windows: usize) -> Self {
        let bins = (0..num_windows)
            .map(|_| Mutex::new(WindowBin { sum_ratio: 0.0, sum_sq: 0.0, count: 0 }))
            .collect();
        WindowStats { bins }
    }

    /// Record a ratio for a specific window position.
    #[inline]
    pub fn record(&self, window_idx: usize, ratio: f64) {
        if ratio.is_finite() {
            let mut bin = self.bins[window_idx].lock().expect("WindowStats bin mutex poisoned");
            bin.sum_ratio += ratio;
            bin.sum_sq += ratio * ratio;
            bin.count += 1;
        }
    }

    /// Get per-window (mean, std_err) for plotting.
    pub fn get_window_data(&self) -> Vec<(f64, f64)> {
        self.bins.iter().map(|b| {
            let bin = b.lock().expect("WindowStats bin mutex poisoned");
            if bin.count == 0 {
                (f64::NAN, f64::NAN)
            } else {
                let mean = bin.sum_ratio / bin.count as f64;
                let se = if bin.count > 1 {
                    let variance = (bin.sum_sq - bin.sum_ratio * bin.sum_ratio / bin.count as f64) / (bin.count - 1) as f64;
                    (variance.max(0.0) / bin.count as f64).sqrt()
                } else {
                    0.0
                };
                (mean, se)
            }
        }).collect()
    }
}

/// Minimal seeded PRNG (SplitMix64) for reproducible bootstrap resampling —
/// avoids pulling in the `rand` crate. Not cryptographic.
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform integer in `[0, bound)` (modulo bias is negligible for the
    /// gene/site counts used in bootstrap resampling).
    #[inline]
    pub fn below(&mut self, bound: usize) -> usize {
        (self.next_u64() % bound as u64) as usize
    }
}

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
fn erf(x: f64) -> f64 {
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

#[cfg(test)]
mod tests {
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
}
