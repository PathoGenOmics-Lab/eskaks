use crate::stats::{SummaryStats, WindowStats};
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufWriter, Write};

mod bars;
mod histogram;
mod window;
#[cfg(test)]
mod tests;

pub use bars::{group_bar_svg, lineage_bar_svg};
pub use histogram::histogram_svg;
pub use window::window_plot_svg;

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

/// XML-escape a user-controlled label before placing it in an SVG `<text>` element.
/// Without this, a genome/lineage/group name containing `&`, `<`, or `>` produces
/// ill-formed XML and the whole `.svg` renders blank in conformant viewers.
fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;")
}

/// Data for group average bar chart.
#[derive(Clone)]
pub struct GroupPlotData {
    pub label: String,
    pub mean: f64,
    pub ci_low: f64,
    pub ci_high: f64,
}

/// Data for lineage bar chart.
pub struct LineagePlotData {
    pub genome: String,
    pub lineage: String,
    pub ratio: f64,
}

