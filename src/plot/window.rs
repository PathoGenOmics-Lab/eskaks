//! Sliding-window dN/dS line plot SVG.

use super::*;

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

