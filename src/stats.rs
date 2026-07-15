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
}
