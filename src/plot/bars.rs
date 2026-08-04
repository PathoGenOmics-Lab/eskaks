//! Group-average and lineage bar-chart SVGs.

use super::*;

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
            cx, MARGIN_TOP + PLOT_H + 20.0, cx, MARGIN_TOP + PLOT_H + 20.0, xml_escape(&group.label));
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
        let label = format!("{}/{}", d.genome, d.lineage);
        let label_short: String = if label.len() > 16 {
            format!("{}..", label.chars().take(14).collect::<String>())
        } else {
            label
        };
        let _ = writeln!(svg,
            r#"<text x="{:.1}" y="{:.1}" class="tick-label" text-anchor="end" transform="rotate(-45,{:.1},{:.1})" font-size="8">{}</text>"#,
            cx, MARGIN_TOP + PLOT_H + 15.0, cx, MARGIN_TOP + PLOT_H + 15.0, xml_escape(&label_short));
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

