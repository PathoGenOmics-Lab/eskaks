//! Manhattan-style SVG plots for the per-gene scan.

use super::*;

/// Generate an SVG Manhattan-style plot of pN/pS per gene along the genome.
pub fn write_pnps_plot(results: &[GenePnPs], prefix: &str, fdr: f64) -> anyhow::Result<String> {
    use std::fmt::Write as FmtWrite;
    use std::fs::File;
    use std::io::{BufWriter, Write};

    // Color constants — passed as format arguments so they never appear
    // literally inside r#"..."# delimiters (which `"#` would terminate).
    const C_GRID: &str = "#e0e0e0";
    const C_POSITIVE: &str = "#d94a4a";
    const C_PURIFYING: &str = "#4a90d9";
    const C_AXIS: &str = "#333333";

    let plot_path = format!("{}_pnps_manhattan.svg", prefix);

    // Filter out genes with no SNPs or NaN/infinite pN/pS
    let plot_data: Vec<&GenePnPs> = results
        .iter()
        .filter(|r| r.total_snps > 0.0 && r.pn_ps.is_finite())
        .collect();

    // Genes with variation but an undefined ratio (pS = 0, i.e. only
    // nonsynonymous SNPs — often the strongest positive-selection candidates)
    // can't be placed on the ratio axis. They stay in the TSV, but say so
    // rather than dropping them from the plot without notice.
    let undefined = results
        .iter()
        .filter(|r| r.total_snps > 0.0 && !r.pn_ps.is_finite())
        .count();
    if undefined > 0 {
        warn!(
            "{} gene(s) with only nonsynonymous variation (pN/pS = ∞) are not shown on the Manhattan plot; they are in the {}_pnps table.",
            undefined, prefix
        );
    }

    if plot_data.is_empty() {
        info!("No valid data points for pN/pS plot");
        return Ok(plot_path);
    }

    let width = 900.0f64;
    let height = 500.0f64;
    let margin_top = 50.0f64;
    let margin_right = 40.0f64;
    let margin_bottom = 80.0f64;
    let margin_left = 80.0f64;
    let plot_w = width - margin_left - margin_right;
    let plot_h = height - margin_top - margin_bottom;

    let max_pos = plot_data.iter().map(|r| r.genome_start).max().unwrap_or(1) as f64;
    let min_pos = plot_data.iter().map(|r| r.genome_start).min().unwrap_or(0) as f64;
    let pos_range = (max_pos - min_pos).max(1.0);

    let max_y = plot_data
        .iter()
        .map(|r| r.pn_ps)
        .fold(f64::NEG_INFINITY, f64::max)
        .max(1.5)
        * 1.1;

    let to_x = |pos: usize| margin_left + ((pos as f64 - min_pos) / pos_range) * plot_w;
    let to_y = |v: f64| margin_top + plot_h * (1.0 - v / max_y);

    let mut svg = String::with_capacity(4096);

    // Header
    svg.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {w} {h}\" width=\"{w}\" height=\"{h}\">",
        w = width, h = height
    );
    let _ = writeln!(svg, "<style>");
    let _ = writeln!(svg, "  text {{ font-family: sans-serif; fill: {c}; }}", c = C_AXIS);
    let _ = writeln!(svg, "  .title {{ font-size: 16px; font-weight: bold; text-anchor: middle; }}");
    let _ = writeln!(svg, "  .axis-label {{ font-size: 12px; text-anchor: middle; }}");
    let _ = writeln!(svg, "  .tick-label {{ font-size: 10px; }}");
    let _ = writeln!(svg, "</style>");
    let _ = writeln!(svg, "<rect width=\"{w}\" height=\"{h}\" fill=\"white\"/>", w = width, h = height);
    let _ = writeln!(
        svg,
        "<text x=\"{cx}\" y=\"30\" class=\"title\">pN/pS per Gene (Manhattan Plot)</text>",
        cx = width / 2.0
    );

    // Y-axis grid lines and labels
    let num_y_ticks = 5;
    for i in 0..=num_y_ticks {
        let frac = i as f64 / num_y_ticks as f64;
        let val = max_y * frac;
        let y = to_y(val);
        let _ = writeln!(
            svg,
            "<line x1=\"{}\" y1=\"{:.1}\" x2=\"{}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"0.5\"/>",
            margin_left, y, margin_left + plot_w, y, C_GRID
        );
        let _ = writeln!(
            svg,
            "<text x=\"{x}\" y=\"{y:.1}\" class=\"tick-label\" text-anchor=\"end\" dominant-baseline=\"middle\">{val:.2}</text>",
            x = margin_left - 8.0, y = y, val = val
        );
    }

    // Neutral line at pN/pS = 1.0
    if 1.0 <= max_y {
        let y1 = to_y(1.0);
        let _ = writeln!(
            svg,
            "<line x1=\"{}\" y1=\"{:.1}\" x2=\"{}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"1\" stroke-dasharray=\"6,3\"/>",
            margin_left, y1, margin_left + plot_w, y1, C_POSITIVE
        );
        let _ = writeln!(
            svg,
            "<text x=\"{x:.1}\" y=\"{y:.1}\" class=\"tick-label\" fill=\"{c}\">pN/pS = 1</text>",
            x = margin_left + plot_w + 3.0, y = y1 + 3.0, c = C_POSITIVE
        );
    }

    // Data points
    for r in &plot_data {
        let x = to_x(r.genome_start);
        let y = to_y(r.pn_ps);
        let color = if r.pn_ps < 1.0 { C_PURIFYING } else { C_POSITIVE };
        let radius = r.total_snps.sqrt().clamp(2.0, 8.0);
        // Outline genes significant in the per-gene neutrality test (BH-FDR).
        let sig = r.q_value.is_finite() && r.q_value < fdr;
        let stroke = if sig {
            " stroke=\"#000000\" stroke-width=\"1.5\""
        } else {
            ""
        };
        let _ = writeln!(
            svg,
            "<circle cx=\"{x:.1}\" cy=\"{y:.1}\" r=\"{r:.1}\" fill=\"{c}\" opacity=\"0.7\"{stroke}>",
            x = x, y = y, r = radius, c = color, stroke = stroke
        );
        let _ = writeln!(
            svg,
            "  <title>{name}: pN/pS={ratio:.4} ({s}S/{n}N SNPs), q={q}</title>",
            name = r.name, ratio = r.pn_ps, s = r.syn_snps, n = r.nonsyn_snps,
            q = format_pval(r.q_value)
        );
        svg.push_str("</circle>\n");
    }

    // Axes
    let _ = writeln!(
        svg,
        "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{}\" stroke-width=\"1.5\"/>",
        margin_left, margin_top, margin_left, margin_top + plot_h, C_AXIS
    );
    let _ = writeln!(
        svg,
        "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{}\" stroke-width=\"1.5\"/>",
        margin_left, margin_top + plot_h, margin_left + plot_w, margin_top + plot_h, C_AXIS
    );

    // Axis labels
    let _ = writeln!(
        svg,
        "<text x=\"{x:.1}\" y=\"{y:.1}\" class=\"axis-label\">Genome Position</text>",
        x = margin_left + plot_w / 2.0, y = height - 10.0
    );
    let _ = writeln!(
        svg,
        "<text x=\"15\" y=\"{y:.1}\" class=\"axis-label\" transform=\"rotate(-90,15,{y:.1})\">pN/pS</text>",
        y = margin_top + plot_h / 2.0
    );

    // Legend
    let _ = writeln!(
        svg,
        "<rect x=\"{x:.1}\" y=\"40\" width=\"12\" height=\"12\" fill=\"{c}\"/>",
        x = width - 200.0, c = C_PURIFYING
    );
    let _ = writeln!(
        svg,
        "<text x=\"{x:.1}\" y=\"50\" class=\"tick-label\">Purifying (pN/pS &lt; 1)</text>",
        x = width - 184.0
    );
    let _ = writeln!(
        svg,
        "<rect x=\"{x:.1}\" y=\"56\" width=\"12\" height=\"12\" fill=\"{c}\"/>",
        x = width - 200.0, c = C_POSITIVE
    );
    let _ = writeln!(
        svg,
        "<text x=\"{x:.1}\" y=\"66\" class=\"tick-label\">Positive (pN/pS &ge; 1)</text>",
        x = width - 184.0
    );

    svg.push_str("</svg>\n");

    let mut file = BufWriter::new(File::create(&plot_path)?);
    file.write_all(svg.as_bytes())?;

    Ok(plot_path)
}

/// Write a −log10(p) Manhattan plot of the per-gene neutrality test, with a
/// Benjamini-Hochberg significance line at the given FDR. Genes significant at
/// that FDR are drawn in red; the rest in grey. Returns the path (an empty file
/// is skipped when no gene has a finite p-value, e.g. under --af-weighted).
pub fn write_pvalue_manhattan(results: &[GenePnPs], prefix: &str, fdr: f64) -> anyhow::Result<String> {
    use std::fmt::Write as FmtWrite;
    use std::fs::File;
    use std::io::{BufWriter, Write};

    const C_GRID: &str = "#e0e0e0";
    const C_SIG: &str = "#d94a4a";
    const C_NS: &str = "#9aa0a6";
    const C_AXIS: &str = "#333333";
    const C_LINE: &str = "#2a7f3f";

    let plot_path = format!("{}_pvalue_manhattan.svg", prefix);
    let data: Vec<&GenePnPs> = results.iter().filter(|r| r.p_value.is_finite()).collect();
    if data.is_empty() {
        info!("No p-values to plot (neutrality test not run)");
        return Ok(plot_path);
    }

    // Benjamini-Hochberg significance threshold p*: the largest p_(i) with
    // p_(i) <= (i/m)·fdr. Points above -log10(p*) are significant.
    let mut sorted: Vec<f64> = data.iter().map(|r| r.p_value).collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = sorted.len();
    let mut p_star = f64::NAN;
    for (i, &p) in sorted.iter().enumerate() {
        if p <= ((i + 1) as f64 / m as f64) * fdr {
            p_star = p;
        }
    }

    let neglog = |p: f64| -(p.max(1e-300).log10());
    let width = 900.0f64;
    let height = 500.0f64;
    let (mt, mr, mb, ml) = (50.0f64, 40.0f64, 80.0f64, 80.0f64);
    let plot_w = width - ml - mr;
    let plot_h = height - mt - mb;

    let max_pos = data.iter().map(|r| r.genome_start).max().unwrap_or(1) as f64;
    let min_pos = data.iter().map(|r| r.genome_start).min().unwrap_or(0) as f64;
    let pos_range = (max_pos - min_pos).max(1.0);
    let max_y = data
        .iter()
        .map(|r| neglog(r.p_value))
        .fold(f64::NEG_INFINITY, f64::max)
        .max(neglog(p_star).max(1.3))
        * 1.1;

    let to_x = |pos: usize| ml + ((pos as f64 - min_pos) / pos_range) * plot_w;
    let to_y = |v: f64| mt + plot_h * (1.0 - v / max_y);

    let mut svg = String::with_capacity(4096);
    svg.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    let _ = writeln!(svg, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {w} {h}\" width=\"{w}\" height=\"{h}\">", w = width, h = height);
    let _ = writeln!(svg, "<style> text {{ font-family: sans-serif; fill: {c}; }} .t {{ font-size:16px; font-weight:bold; text-anchor:middle; }} .a {{ font-size:12px; text-anchor:middle; }} .l {{ font-size:10px; }}</style>", c = C_AXIS);
    let _ = writeln!(svg, "<rect width=\"{w}\" height=\"{h}\" fill=\"white\"/>", w = width, h = height);
    let _ = writeln!(svg, "<text x=\"{cx}\" y=\"30\" class=\"t\">Per-gene neutrality test (−log10 p)</text>", cx = width / 2.0);

    // Y grid + labels
    for i in 0..=5 {
        let val = max_y * i as f64 / 5.0;
        let y = to_y(val);
        let _ = writeln!(svg, "<line x1=\"{}\" y1=\"{:.1}\" x2=\"{}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"0.5\"/>", ml, y, ml + plot_w, y, C_GRID);
        let _ = writeln!(svg, "<text x=\"{x}\" y=\"{y:.1}\" class=\"l\" text-anchor=\"end\" dominant-baseline=\"middle\">{val:.1}</text>", x = ml - 8.0, y = y, val = val);
    }

    // BH significance line
    if p_star.is_finite() {
        let y = to_y(neglog(p_star));
        let _ = writeln!(svg, "<line x1=\"{}\" y1=\"{:.1}\" x2=\"{}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"1\" stroke-dasharray=\"6,3\"/>", ml, y, ml + plot_w, y, C_LINE);
        let _ = writeln!(svg, "<text x=\"{x:.1}\" y=\"{y:.1}\" class=\"l\" fill=\"{c}\">BH FDR {fdr}</text>", x = ml + plot_w + 3.0, y = y + 3.0, c = C_LINE, fdr = fdr);
    }

    // Points
    for r in &data {
        let x = to_x(r.genome_start);
        let y = to_y(neglog(r.p_value));
        let sig = r.q_value.is_finite() && r.q_value < fdr;
        let color = if sig { C_SIG } else { C_NS };
        let radius = r.total_snps.sqrt().clamp(2.0, 8.0);
        let _ = writeln!(svg, "<circle cx=\"{x:.1}\" cy=\"{y:.1}\" r=\"{r:.1}\" fill=\"{c}\" opacity=\"0.75\">", x = x, y = y, r = radius, c = color);
        let _ = writeln!(svg, "  <title>{name}: p={p}, q={q}, pN/pS={ratio:.3}</title>", name = r.name, p = format_pval(r.p_value), q = format_pval(r.q_value), ratio = r.pn_ps);
        svg.push_str("</circle>\n");
    }

    // Axes + labels
    let _ = writeln!(svg, "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{}\" stroke-width=\"1.5\"/>", ml, mt, ml, mt + plot_h, C_AXIS);
    let _ = writeln!(svg, "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{}\" stroke-width=\"1.5\"/>", ml, mt + plot_h, ml + plot_w, mt + plot_h, C_AXIS);
    let _ = writeln!(svg, "<text x=\"{x:.1}\" y=\"{y:.1}\" class=\"a\">Genome Position</text>", x = ml + plot_w / 2.0, y = height - 10.0);
    let _ = writeln!(svg, "<text x=\"15\" y=\"{y:.1}\" class=\"a\" transform=\"rotate(-90,15,{y:.1})\">−log10(p)</text>", y = mt + plot_h / 2.0);

    svg.push_str("</svg>\n");
    let mut file = BufWriter::new(File::create(&plot_path)?);
    file.write_all(svg.as_bytes())?;
    Ok(plot_path)
}

