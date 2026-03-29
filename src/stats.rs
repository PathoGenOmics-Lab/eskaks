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
        let mut shared = self.floats.lock().unwrap();
        shared.merge(local);
    }

    /// Print summary to stderr.
    pub fn print_summary(&self) {
        let total = self.total_count.load(Ordering::Relaxed);
        let nans = self.nan_count.load(Ordering::Relaxed);
        let valid = total.saturating_sub(nans);
        let floats = self.floats.lock().unwrap();

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
        if floats.finite_ratio_count > 0 {
            let mean_ratio = floats.sum_ratio / floats.finite_ratio_count as f64;
            eprintln!("  dN/dS: min={:.6}  max={:.6}  mean={:.6}", floats.min_ratio, floats.max_ratio, mean_ratio);
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
            let mut bin = self.bins[window_idx].lock().unwrap();
            bin.sum_ratio += ratio;
            bin.sum_sq += ratio * ratio;
            bin.count += 1;
        }
    }

    /// Get per-window (mean, std_err) for plotting.
    pub fn get_window_data(&self) -> Vec<(f64, f64)> {
        self.bins.iter().map(|b| {
            let bin = b.lock().unwrap();
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
