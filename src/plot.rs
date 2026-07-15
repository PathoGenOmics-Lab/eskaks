use crate::stats::{SummaryStats, WindowStats};
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufWriter, Write};

// SVG canvas dimensions
const WIDTH: f64 = 800.0;
const HEIGHT: f64 = 500.0;
const MARGIN_TOP: f64 = 50.0;
const MARGIN_RIGHT: f64 = 40.0;
const MARGIN_BOTTOM: f64 = 80.0;
const MARGIN_LEFT: f64 = 80.0;

const PLOT_W: f64 = WIDTH - MARGIN_LEFT - MARGIN_RIGHT;
const PLOT_H: f64 = HEIGHT - MARGIN_TOP - MARGIN_BOTTOM;

const COLOR_PURIFYING: &str = "#4a90d9";
const COLOR_POSITIVE: &str = "#d94a4a";
const COLOR_LINE: &str = "#2c3e50";
const COLOR_CI_BAND: &str = "rgba(74,144,217,0.2)";
const COLOR_GRID: &str = "#e0e0e0";
const COLOR_AXIS: &str = "#333333";

fn svg_header(title: &str) -> String {
    format!(
        r#"<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">
<style>
  text {{ font-family: sans-serif; fill: {axis}; }}
  .title {{ font-size: 16px; font-weight: bold; text-anchor: middle; }}
  .axis-label {{ font-size: 12px; text-anchor: middle; }}
  .tick-label {{ font-size: 10px; }}
  .grid {{ stroke: {grid}; stroke-width: 0.5; }}
</style>
<rect width="{w}" height="{h}" fill="white"/>
<text x="{cx}" y="30" class="title">{title}</text>
"#,
        w = WIDTH, h = HEIGHT, cx = WIDTH / 2.0,
        axis = COLOR_AXIS, grid = COLOR_GRID, title = title
    )
}

fn svg_footer() -> &'static str {
    "</svg>\n"
}

/// Write SVG content to a file.
fn write_svg(path: &str, content: &str) -> std::io::Result<()> {
    let mut file = BufWriter::new(
        File::create(path)
            .map_err(|e| std::io::Error::new(e.kind(), format!("Cannot create '{}': {}", path, e)))?
    );
    file.write_all(content.as_bytes())?;
    Ok(())
}

/// Generate a histogram SVG of dN/dS distribution.
pub fn histogram_svg(stats: &SummaryStats, path: &str) -> std::io::Result<()> {
    let bins = stats.get_histogram();
    let max_count = bins.iter().map(|(_, c)| *c).max().unwrap_or(1).max(1);
    let num_bins = bins.len();
    let bar_width = PLOT_W / num_bins as f64 * 0.8;
    let gap = PLOT_W / num_bins as f64 * 0.2;

    let mut svg = svg_header("dN/dS Ratio Distribution");

    // Y-axis grid lines and labels — use integer steps to avoid duplicate "0" ticks
    let num_y_ticks = max_count.clamp(1, 5);
    let step = max_count.div_ceil(num_y_ticks).max(1);
    let mut tick_val = 0usize;
    while tick_val <= max_count {
        let frac = tick_val as f64 / max_count as f64;
        let y = MARGIN_TOP + PLOT_H * (1.0 - frac);
        let _ = writeln!(svg,
            r#"<line x1="{}" y1="{:.1}" x2="{}" y2="{:.1}" class="grid"/>"#,
            MARGIN_LEFT, y, MARGIN_LEFT + PLOT_W, y);
        let _ = writeln!(svg,
            r#"<text x="{}" y="{:.1}" class="tick-label" text-anchor="end" dominant-baseline="middle">{}</text>"#,
            MARGIN_LEFT - 8.0, y, tick_val);
        tick_val += step;
    }

    // Bars
    for (i, (label, count)) in bins.iter().enumerate() {
        let bar_height = (*count as f64 / max_count as f64) * PLOT_H;
        let x = MARGIN_LEFT + i as f64 * (bar_width + gap) + gap / 2.0;
        let y = MARGIN_TOP + PLOT_H - bar_height;
        let color = if i < 5 { COLOR_PURIFYING } else { COLOR_POSITIVE };

        let _ = writeln!(svg,
            r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="{}" rx="2"/>"#,
            x, y, bar_width, bar_height, color);

        // Bar label (count)
        if *count > 0 {
            let _ = writeln!(svg,
                r#"<text x="{:.1}" y="{:.1}" class="tick-label" text-anchor="middle">{}</text>"#,
                x + bar_width / 2.0, y - 5.0, count);
        }

        // X-axis label
        let _ = writeln!(svg,
            r#"<text x="{:.1}" y="{:.1}" class="tick-label" text-anchor="middle" transform="rotate(-30,{:.1},{:.1})">{}</text>"#,
            x + bar_width / 2.0, MARGIN_TOP + PLOT_H + 20.0,
            x + bar_width / 2.0, MARGIN_TOP + PLOT_H + 20.0,
            label);
    }

    // Axes
    let _ = writeln!(svg,
        r#"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="1.5"/>"#,
        MARGIN_LEFT, MARGIN_TOP, MARGIN_LEFT, MARGIN_TOP + PLOT_H, COLOR_AXIS);
    let _ = writeln!(svg,
        r#"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="1.5"/>"#,
        MARGIN_LEFT, MARGIN_TOP + PLOT_H, MARGIN_LEFT + PLOT_W, MARGIN_TOP + PLOT_H, COLOR_AXIS);

    // Axis labels
    let _ = writeln!(svg,
        r#"<text x="{:.1}" y="{:.1}" class="axis-label">dN/dS Ratio</text>"#,
        MARGIN_LEFT + PLOT_W / 2.0, HEIGHT - 10.0);
    let _ = writeln!(svg,
        r#"<text x="15" y="{:.1}" class="axis-label" transform="rotate(-90,15,{:.1})">Count</text>"#,
        MARGIN_TOP + PLOT_H / 2.0, MARGIN_TOP + PLOT_H / 2.0);

    // Legend
    let _ = writeln!(svg,
        r#"<rect x="{:.1}" y="40" width="12" height="12" fill="{}"/>"#,
        WIDTH - 180.0, COLOR_PURIFYING);
    let _ = writeln!(svg,
        r#"<text x="{:.1}" y="50" class="tick-label">Purifying (dN/dS &lt; 1)</text>"#,
        WIDTH - 164.0);
    let _ = writeln!(svg,
        r#"<rect x="{:.1}" y="56" width="12" height="12" fill="{}"/>"#,
        WIDTH - 180.0, COLOR_POSITIVE);
    let _ = writeln!(svg,
        r#"<text x="{:.1}" y="66" class="tick-label">Positive (dN/dS &ge; 1)</text>"#,
        WIDTH - 164.0);

    svg.push_str(svg_footer());
    write_svg(path, &svg)
}

/// Generate a sliding window line plot SVG.
pub fn window_plot_svg(window_stats: &WindowStats, path: &str, window_size: usize, window_step: usize) -> std::io::Result<()> {
    let data = window_stats.get_window_data();
    if data.is_empty() {
        return Ok(());
    }

    let mut svg = svg_header("Sliding Window dN/dS");

    // Filter valid data points
    let valid_data: Vec<(usize, f64, f64)> = data.iter().enumerate()
        .filter(|(_, (m, _))| m.is_finite())
        .map(|(i, (m, se))| (i, *m, *se))
        .collect();

    if valid_data.is_empty() {
        svg.push_str(svg_footer());
        return write_svg(path, &svg);
    }

    // Compute Y range
    let max_y = valid_data.iter().map(|(_, m, se)| m + 1.96 * se)
        .fold(f64::NEG_INFINITY, f64::max)
        .max(1.5);
    let min_y = 0.0_f64;
    let y_range = max_y - min_y;

    let num_windows = data.len();
    let x_scale = PLOT_W / (num_windows.max(1) - 1).max(1) as f64;

    let to_x = |i: usize| MARGIN_LEFT + i as f64 * x_scale;
    let to_y = |v: f64| MARGIN_TOP + PLOT_H * (1.0 - (v - min_y) / y_range);

    // Y-axis grid
    let num_y_ticks = 5;
    for i in 0..=num_y_ticks {
        let frac = i as f64 / num_y_ticks as f64;
        let val = min_y + y_range * frac;
        let y = to_y(val);
        let _ = writeln!(svg,
            r#"<line x1="{}" y1="{:.1}" x2="{}" y2="{:.1}" class="grid"/>"#,
            MARGIN_LEFT, y, MARGIN_LEFT + PLOT_W, y);
        let _ = writeln!(svg,
            r#"<text x="{}" y="{:.1}" class="tick-label" text-anchor="end" dominant-baseline="middle">{:.2}</text>"#,
            MARGIN_LEFT - 8.0, y, val);
    }

    // Neutral line at dN/dS = 1.0
    if 1.0 >= min_y && 1.0 <= max_y {
        let y1 = to_y(1.0);
        let _ = writeln!(svg,
            r#"<line x1="{}" y1="{:.1}" x2="{}" y2="{:.1}" stroke="{}" stroke-width="1" stroke-dasharray="6,3"/>"#,
            MARGIN_LEFT, y1, MARGIN_LEFT + PLOT_W, y1, COLOR_POSITIVE);
        let _ = writeln!(svg,
            r#"<text x="{:.1}" y="{:.1}" class="tick-label" fill="{}">dN/dS = 1</text>"#,
            MARGIN_LEFT + PLOT_W + 3.0, y1 + 3.0, COLOR_POSITIVE);
    }

    // CI band (shaded area)
    if valid_data.len() > 1 {
        let mut path_d = String::new();
        // Upper bound forward
        for (i, (idx, m, se)) in valid_data.iter().enumerate() {
            let x = to_x(*idx);
            let y = to_y((m + 1.96 * se).min(max_y));
            if i == 0 { let _ = write!(path_d, "M{:.1},{:.1}", x, y); }
            else { let _ = write!(path_d, " L{:.1},{:.1}", x, y); }
        }
        // Lower bound backward
        for (idx, m, se) in valid_data.iter().rev() {
            let x = to_x(*idx);
            let y = to_y((m - 1.96 * se).max(min_y));
            let _ = write!(path_d, " L{:.1},{:.1}", x, y);
        }
        let _ = write!(path_d, " Z");
        let _ = writeln!(svg, r#"<path d="{}" fill="{}" stroke="none"/>"#, path_d, COLOR_CI_BAND);
    }

    // Mean line
    if valid_data.len() > 1 {
        let mut path_d = String::new();
        for (i, (idx, m, _)) in valid_data.iter().enumerate() {
            let x = to_x(*idx);
            let y = to_y(*m);
            if i == 0 { let _ = write!(path_d, "M{:.1},{:.1}", x, y); }
            else { let _ = write!(path_d, " L{:.1},{:.1}", x, y); }
        }
        let _ = writeln!(svg,
            r#"<path d="{}" fill="none" stroke="{}" stroke-width="2"/>"#,
            path_d, COLOR_LINE);
    }

    // Axes
    let _ = writeln!(svg,
        r#"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="1.5"/>"#,
        MARGIN_LEFT, MARGIN_TOP, MARGIN_LEFT, MARGIN_TOP + PLOT_H, COLOR_AXIS);
    let _ = writeln!(svg,
        r#"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="1.5"/>"#,
        MARGIN_LEFT, MARGIN_TOP + PLOT_H, MARGIN_LEFT + PLOT_W, MARGIN_TOP + PLOT_H, COLOR_AXIS);

    // X-axis labels (every ~5 ticks)
    let tick_interval = (num_windows / 5).max(1);
    for i in (0..num_windows).step_by(tick_interval) {
        let x = to_x(i);
        let codon_pos = i * window_step + 1;
        let _ = writeln!(svg,
            r#"<text x="{:.1}" y="{:.1}" class="tick-label" text-anchor="middle">{}</text>"#,
            x, MARGIN_TOP + PLOT_H + 18.0, codon_pos);
    }

    // Axis labels
    let _ = writeln!(svg,
        r#"<text x="{:.1}" y="{:.1}" class="axis-label">Codon Position (window size={})</text>"#,
        MARGIN_LEFT + PLOT_W / 2.0, HEIGHT - 10.0, window_size);
    let _ = writeln!(svg,
        r#"<text x="15" y="{:.1}" class="axis-label" transform="rotate(-90,15,{:.1})">Mean dN/dS</text>"#,
        MARGIN_TOP + PLOT_H / 2.0, MARGIN_TOP + PLOT_H / 2.0);

    svg.push_str(svg_footer());
    write_svg(path, &svg)
}

/// Data for group average bar chart.
#[derive(Clone)]
pub struct GroupPlotData {
    pub label: String,
    pub mean: f64,
    pub ci_low: f64,
    pub ci_high: f64,
}

/// Generate a bar chart with error bars for group average dN/dS.
pub fn group_bar_svg(groups: &[GroupPlotData], path: &str) -> std::io::Result<()> {
    if groups.is_empty() {
        return Ok(());
    }

    let mut svg = svg_header("Group Average dN/dS");

    let num_bars = groups.len();
    let bar_width = (PLOT_W / num_bars as f64 * 0.7).min(60.0);
    let bar_spacing = PLOT_W / num_bars as f64;

    // Filter non-finite values to avoid NaN/inf poisoning the scale
    let max_y = groups.iter()
        .map(|g| {
            if g.ci_high.is_finite() { g.ci_high }
            else if g.mean.is_finite() { g.mean }
            else { 0.0 }
        })
        .fold(f64::NEG_INFINITY, f64::max)
        .max(0.1) * 1.2;

    let to_y = |v: f64| MARGIN_TOP + PLOT_H * (1.0 - v / max_y);

    // Y-axis grid
    let num_y_ticks = 5;
    for i in 0..=num_y_ticks {
        let frac = i as f64 / num_y_ticks as f64;
        let val = max_y * frac;
        let y = to_y(val);
        let _ = writeln!(svg,
            r#"<line x1="{}" y1="{:.1}" x2="{}" y2="{:.1}" class="grid"/>"#,
            MARGIN_LEFT, y, MARGIN_LEFT + PLOT_W, y);
        let _ = writeln!(svg,
            r#"<text x="{}" y="{:.1}" class="tick-label" text-anchor="end" dominant-baseline="middle">{:.3}</text>"#,
            MARGIN_LEFT - 8.0, y, val);
    }

    // Neutral line at 1.0
    if 1.0 <= max_y {
        let y1 = to_y(1.0);
        let _ = writeln!(svg,
            r#"<line x1="{}" y1="{:.1}" x2="{}" y2="{:.1}" stroke="{}" stroke-width="1" stroke-dasharray="6,3"/>"#,
            MARGIN_LEFT, y1, MARGIN_LEFT + PLOT_W, y1, COLOR_POSITIVE);
    }

    // Bars and error bars
    for (i, group) in groups.iter().enumerate() {
        let cx = MARGIN_LEFT + (i as f64 + 0.5) * bar_spacing;
        let x = cx - bar_width / 2.0;
        let mean_clamped = if group.mean.is_finite() { group.mean.max(0.0) } else { max_y };
        let bar_h = PLOT_H * (mean_clamped / max_y);
        let y = to_y(mean_clamped);
        let color = if group.mean.is_finite() && group.mean < 1.0 { COLOR_PURIFYING } else { COLOR_POSITIVE };

        // Bar
        let _ = writeln!(svg,
            r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="{}" rx="2"/>"#,
            x, y, bar_width, bar_h, color);

        // Label for infinite values
        if !group.mean.is_finite() {
            let _ = writeln!(svg,
                r#"<text x="{:.1}" y="{:.1}" class="tick-label" text-anchor="middle" font-weight="bold">inf</text>"#,
                cx, y - 5.0);
        }

        // Error bars (whiskers)
        if group.ci_low.is_finite() && group.ci_high.is_finite() {
            let y_low = to_y(group.ci_low.max(0.0));
            let y_high = to_y(group.ci_high.min(max_y));
            let whisker_w = bar_width * 0.3;
            // Vertical line
            let _ = writeln!(svg,
                r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="1.5"/>"#,
                cx, y_high, cx, y_low, COLOR_AXIS);
            // Top whisker
            let _ = writeln!(svg,
                r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="1.5"/>"#,
                cx - whisker_w / 2.0, y_high, cx + whisker_w / 2.0, y_high, COLOR_AXIS);
            // Bottom whisker
            let _ = writeln!(svg,
                r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="1.5"/>"#,
                cx - whisker_w / 2.0, y_low, cx + whisker_w / 2.0, y_low, COLOR_AXIS);
        }

        // X-axis label
        let _ = writeln!(svg,
            r#"<text x="{:.1}" y="{:.1}" class="tick-label" text-anchor="middle" transform="rotate(-30,{:.1},{:.1})">{}</text>"#,
            cx, MARGIN_TOP + PLOT_H + 20.0, cx, MARGIN_TOP + PLOT_H + 20.0, group.label);
    }

    // Axes
    let _ = writeln!(svg,
        r#"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="1.5"/>"#,
        MARGIN_LEFT, MARGIN_TOP, MARGIN_LEFT, MARGIN_TOP + PLOT_H, COLOR_AXIS);
    let _ = writeln!(svg,
        r#"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="1.5"/>"#,
        MARGIN_LEFT, MARGIN_TOP + PLOT_H, MARGIN_LEFT + PLOT_W, MARGIN_TOP + PLOT_H, COLOR_AXIS);

    // Axis labels
    let _ = writeln!(svg,
        r#"<text x="{:.1}" y="{:.1}" class="axis-label">Group Pair</text>"#,
        MARGIN_LEFT + PLOT_W / 2.0, HEIGHT - 10.0);
    let _ = writeln!(svg,
        r#"<text x="15" y="{:.1}" class="axis-label" transform="rotate(-90,15,{:.1})">Mean dN/dS</text>"#,
        MARGIN_TOP + PLOT_H / 2.0, MARGIN_TOP + PLOT_H / 2.0);

    svg.push_str(svg_footer());
    write_svg(path, &svg)
}

/// Data for lineage bar chart.
pub struct LineagePlotData {
    pub genome: String,
    pub lineage: String,
    pub ratio: f64,
}

/// Generate a bar chart for lineage dN/dS ratios.
pub fn lineage_bar_svg(data: &[LineagePlotData], path: &str) -> std::io::Result<()> {
    if data.is_empty() {
        return Ok(());
    }

    let mut svg = svg_header("Lineage dN/dS Ratios");

    let num_bars = data.len().min(50); // limit bars for readability
    let bar_width = (PLOT_W / num_bars as f64 * 0.8).min(40.0);
    let bar_spacing = PLOT_W / num_bars as f64;

    let max_y = data.iter().take(num_bars)
        .map(|d| if d.ratio.is_finite() { d.ratio } else { 0.0 })
        .fold(f64::NEG_INFINITY, f64::max)
        .max(0.1) * 1.2;

    let to_y = |v: f64| MARGIN_TOP + PLOT_H * (1.0 - v / max_y);

    // Y-axis grid
    let num_y_ticks = 5;
    for i in 0..=num_y_ticks {
        let frac = i as f64 / num_y_ticks as f64;
        let val = max_y * frac;
        let y = to_y(val);
        let _ = writeln!(svg,
            r#"<line x1="{}" y1="{:.1}" x2="{}" y2="{:.1}" class="grid"/>"#,
            MARGIN_LEFT, y, MARGIN_LEFT + PLOT_W, y);
        let _ = writeln!(svg,
            r#"<text x="{}" y="{:.1}" class="tick-label" text-anchor="end" dominant-baseline="middle">{:.2}</text>"#,
            MARGIN_LEFT - 8.0, y, val);
    }

    // Neutral line
    if 1.0 <= max_y {
        let y1 = to_y(1.0);
        let _ = writeln!(svg,
            r#"<line x1="{}" y1="{:.1}" x2="{}" y2="{:.1}" stroke="{}" stroke-width="1" stroke-dasharray="6,3"/>"#,
            MARGIN_LEFT, y1, MARGIN_LEFT + PLOT_W, y1, COLOR_POSITIVE);
    }

    for (i, d) in data.iter().take(num_bars).enumerate() {
        let cx = MARGIN_LEFT + (i as f64 + 0.5) * bar_spacing;
        let x = cx - bar_width / 2.0;
        let ratio = if d.ratio.is_finite() { d.ratio.max(0.0) } else { 0.0 };
        let bar_h = PLOT_H * (ratio / max_y);
        let y = to_y(ratio);
        let color = if ratio < 1.0 { COLOR_PURIFYING } else { COLOR_POSITIVE };

        let _ = writeln!(svg,
            r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="{}" rx="2"/>"#,
            x, y, bar_width, bar_h, color);

        // Label (truncated for readability)
        let label = format!("{}/{}", &d.genome, &d.lineage);
        let label_short: String = if label.len() > 16 {
            format!("{}..", label.chars().take(14).collect::<String>())
        } else {
            label
        };
        let _ = writeln!(svg,
            r#"<text x="{:.1}" y="{:.1}" class="tick-label" text-anchor="end" transform="rotate(-45,{:.1},{:.1})" font-size="8">{}</text>"#,
            cx, MARGIN_TOP + PLOT_H + 15.0, cx, MARGIN_TOP + PLOT_H + 15.0, label_short);
    }

    // Axes
    let _ = writeln!(svg,
        r#"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="1.5"/>"#,
        MARGIN_LEFT, MARGIN_TOP, MARGIN_LEFT, MARGIN_TOP + PLOT_H, COLOR_AXIS);
    let _ = writeln!(svg,
        r#"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="1.5"/>"#,
        MARGIN_LEFT, MARGIN_TOP + PLOT_H, MARGIN_LEFT + PLOT_W, MARGIN_TOP + PLOT_H, COLOR_AXIS);

    let _ = writeln!(svg,
        r#"<text x="{:.1}" y="{:.1}" class="axis-label">Genome / Lineage</text>"#,
        MARGIN_LEFT + PLOT_W / 2.0, HEIGHT - 5.0);
    let _ = writeln!(svg,
        r#"<text x="15" y="{:.1}" class="axis-label" transform="rotate(-90,15,{:.1})">dN/dS Ratio</text>"#,
        MARGIN_TOP + PLOT_H / 2.0, MARGIN_TOP + PLOT_H / 2.0);

    svg.push_str(svg_footer());
    write_svg(path, &svg)
}
