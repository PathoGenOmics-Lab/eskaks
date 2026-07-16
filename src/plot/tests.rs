use super::*;

#[test]
fn xml_escape_neutralizes_metacharacters() {
    assert_eq!(xml_escape("L2&L4"), "L2&amp;L4");
    assert_eq!(xml_escape("a<b>c"), "a&lt;b&gt;c");
    assert_eq!(xml_escape("plain"), "plain");
    // & is escaped first, so an existing entity is not double-mangled into a valid one.
    assert_eq!(xml_escape("x & <y>"), "x &amp; &lt;y&gt;");
}

/// Write to a temp path and return the resulting SVG string.
fn render(f: impl FnOnce(&str) -> std::io::Result<()>) -> String {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("plot.svg");
    let path_str = path.to_str().unwrap();
    f(path_str).unwrap();
    std::fs::read_to_string(&path).unwrap_or_default()
}

/// Structural + numeric-safety invariants every non-empty SVG we emit must hold.
///
/// `NaN` and `inf` never belong in a coordinate/length attribute: a single
/// poisoned attribute makes conformant SVG viewers drop the whole element (or
/// the whole document). The one legitimate "inf" in our output is the group
/// bar's `>inf</text>` overflow label — that is element *text*, never an
/// `attr="inf"`, so the attribute-scoped check below leaves it alone.
fn assert_well_formed(svg: &str) {
    assert!(svg.starts_with("<?xml"), "missing XML prolog: {:?}", &svg[..svg.len().min(40)]);
    assert!(svg.contains("<svg "), "missing <svg> root");
    assert!(svg.trim_end().ends_with("</svg>"), "not closed with </svg>");
    assert!(!svg.contains("NaN"), "NaN leaked into SVG output");
    // Attribute-scoped non-finite check (leaves the intentional >inf</text> label).
    for bad in ["=\"inf\"", "=\"-inf\"", "=\"NaN\""] {
        assert!(!svg.contains(bad), "non-finite value in attribute ({bad})");
    }
}

fn stats_with(ratios: &[f64]) -> SummaryStats {
    let s = SummaryStats::new();
    for &r in ratios {
        // dn/ds are only used for the mean summary, not the histogram bins.
        s.record_pair_atomic(r * 0.5, 0.5, r);
    }
    s
}

// ---- histogram_svg ----------------------------------------------------

#[test]
fn histogram_typical_distribution_is_well_formed() {
    let svg = render(|p| histogram_svg(&stats_with(&[0.1, 0.3, 0.7, 0.95, 1.0, 1.5, 3.0]), p));
    assert_well_formed(&svg);
    assert!(svg.contains("dN/dS Ratio Distribution"));
    assert!(svg.contains("<rect"));
}

#[test]
fn histogram_empty_stats_still_valid() {
    // No pairs recorded: every bin is 0, max_count falls back to 1.
    let svg = render(|p| histogram_svg(&SummaryStats::new(), p));
    assert_well_formed(&svg);
}

#[test]
fn histogram_single_bin_populated() {
    let svg = render(|p| histogram_svg(&stats_with(&[0.5, 0.5, 0.5]), p));
    assert_well_formed(&svg);
}

// ---- window_plot_svg --------------------------------------------------

fn windows_with(rows: &[&[f64]]) -> WindowStats {
    let w = WindowStats::new(rows.len());
    for (i, row) in rows.iter().enumerate() {
        for &r in *row {
            w.record(i, r);
        }
    }
    w
}

#[test]
fn window_typical_series_is_well_formed() {
    let w = windows_with(&[&[0.4, 0.5, 0.6], &[0.9, 1.0, 1.1], &[1.5, 1.6, 1.7], &[0.2, 0.3]]);
    let svg = render(|p| window_plot_svg(&w, p, 30, 3));
    assert_well_formed(&svg);
    assert!(svg.contains("Sliding Window dN/dS"));
    assert!(svg.contains("<path")); // CI band + mean line for >1 valid point
}

#[test]
fn window_empty_writes_nothing() {
    // get_window_data().is_empty() => early Ok(()) without touching the file.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("w.svg");
    window_plot_svg(&WindowStats::new(0), path.to_str().unwrap(), 30, 3).unwrap();
    assert!(!path.exists(), "empty window stats must not create a file");
}

#[test]
fn window_all_empty_bins_is_well_formed() {
    // Windows exist but no ratios recorded: every (mean,se) is (NaN,NaN) and
    // gets filtered out. The header/footer path must still be NaN-free.
    let svg = render(|p| window_plot_svg(&WindowStats::new(4), p, 30, 3));
    assert_well_formed(&svg);
}

#[test]
fn window_single_valid_point_no_line() {
    // Only one window has data: no CI band / mean line (needs >1), but valid SVG.
    let w = windows_with(&[&[], &[1.2, 1.3], &[]]);
    let svg = render(|p| window_plot_svg(&w, p, 30, 3));
    assert_well_formed(&svg);
}

#[test]
fn window_zero_variance_gives_zero_se_not_nan() {
    // A single sample per window => se branch returns 0.0 (never 0/0 = NaN).
    let w = windows_with(&[&[1.0], &[1.0], &[1.0]]);
    let svg = render(|p| window_plot_svg(&w, p, 30, 3));
    assert_well_formed(&svg);
}

// ---- group_bar_svg ----------------------------------------------------

#[test]
fn group_bar_typical_is_well_formed() {
    let groups = vec![
        GroupPlotData { label: "L1 vs L2".into(), mean: 0.4, ci_low: 0.3, ci_high: 0.5 },
        GroupPlotData { label: "L2 vs L4".into(), mean: 1.3, ci_low: 1.0, ci_high: 1.7 },
    ];
    let svg = render(|p| group_bar_svg(&groups, p));
    assert_well_formed(&svg);
    assert!(svg.contains("Group Average dN/dS"));
}

#[test]
fn group_bar_empty_writes_nothing() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("g.svg");
    group_bar_svg(&[], path.to_str().unwrap()).unwrap();
    assert!(!path.exists(), "empty groups must not create a file");
}

#[test]
fn group_bar_infinite_mean_renders_inf_label_not_nan() {
    // pN/pS = ∞ (zero synonymous) must degrade to an "inf" text label and must
    // not poison any coordinate attribute with NaN/inf.
    let groups = vec![
        GroupPlotData { label: "zeroS".into(), mean: f64::INFINITY, ci_low: f64::INFINITY, ci_high: f64::INFINITY },
        GroupPlotData { label: "ok".into(), mean: 0.6, ci_low: 0.4, ci_high: 0.8 },
    ];
    let svg = render(|p| group_bar_svg(&groups, p));
    assert_well_formed(&svg);
    assert!(svg.contains(">inf</text>"), "expected an 'inf' overflow label");
}

#[test]
fn group_bar_escapes_label_metacharacters() {
    let groups = vec![
        GroupPlotData { label: "A&B<C>".into(), mean: 0.5, ci_low: 0.4, ci_high: 0.6 },
    ];
    let svg = render(|p| group_bar_svg(&groups, p));
    assert_well_formed(&svg);
    assert!(svg.contains("A&amp;B&lt;C&gt;"), "label not XML-escaped");
    assert!(!svg.contains("A&B<C>"), "raw metacharacters leaked into SVG");
}

#[test]
fn group_bar_negative_mean_clamps_to_zero() {
    let groups = vec![
        GroupPlotData { label: "neg".into(), mean: -0.5, ci_low: -1.0, ci_high: 0.2 },
    ];
    let svg = render(|p| group_bar_svg(&groups, p));
    assert_well_formed(&svg);
}

// ---- lineage_bar_svg --------------------------------------------------

#[test]
fn lineage_bar_typical_is_well_formed() {
    let data = vec![
        LineagePlotData { genome: "H37Rv".into(), lineage: "L2".into(), ratio: 0.4 },
        LineagePlotData { genome: "H37Rv".into(), lineage: "L4".into(), ratio: 1.6 },
    ];
    let svg = render(|p| lineage_bar_svg(&data, p));
    assert_well_formed(&svg);
    assert!(svg.contains("Lineage dN/dS Ratios"));
}

#[test]
fn lineage_bar_empty_writes_nothing() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("l.svg");
    lineage_bar_svg(&[], path.to_str().unwrap()).unwrap();
    assert!(!path.exists(), "empty lineage data must not create a file");
}

#[test]
fn lineage_bar_infinite_ratio_stays_finite_in_svg() {
    let data = vec![
        LineagePlotData { genome: "g".into(), lineage: "inf".into(), ratio: f64::INFINITY },
        LineagePlotData { genome: "g".into(), lineage: "ok".into(), ratio: 0.7 },
    ];
    let svg = render(|p| lineage_bar_svg(&data, p));
    assert_well_formed(&svg);
}

#[test]
fn lineage_bar_escapes_and_truncates_long_label() {
    // Long label (>16 chars) is truncated with ".." and must still be escaped.
    let data = vec![
        LineagePlotData { genome: "genome_with_a_really_long_name".into(), lineage: "L2&<x>".into(), ratio: 0.9 },
    ];
    let svg = render(|p| lineage_bar_svg(&data, p));
    assert_well_formed(&svg);
    assert!(!svg.contains("L2&<x>"), "raw metacharacters leaked");
}

#[test]
fn lineage_bar_caps_at_fifty_bars() {
    // 60 entries provided; only the first 50 are drawn. Count <rect> data bars
    // (exclude the header white background rect).
    let data: Vec<LineagePlotData> = (0..60)
        .map(|i| LineagePlotData { genome: format!("g{i}"), lineage: "L".into(), ratio: 0.5 })
        .collect();
    let svg = render(|p| lineage_bar_svg(&data, p));
    assert_well_formed(&svg);
    let bar_count = svg.matches(r#"rx="2""#).count(); // data bars carry rx="2"
    assert_eq!(bar_count, 50, "should cap at 50 bars, got {bar_count}");
}
