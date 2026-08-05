//! dN/dS histogram SVG.

use super::*;

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
    let step = ((max_count + num_y_ticks - 1) / num_y_ticks).max(1);
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

