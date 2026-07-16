use crate::codon::extract_group_key;
use crate::models::{DsDn, Model, Z_95_CONFIDENCE};
use crate::plot::GroupPlotData;
use crate::stats::{FloatAccum, SummaryStats, WindowStats};
use crossbeam::channel::unbounded;
use indicatif::{ProgressBar, ProgressStyle};
use log::info;
use rayon::prelude::*;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

mod grouped;
mod pairwise;
mod windows;
#[cfg(test)]
mod tests;

// Public writers (re-exported so callers keep using `output::…`).
pub use grouped::{write_group_average, write_lineage};
pub use pairwise::{write_pairwise, write_pairwise_bootstrap, write_pairwise_tests};
pub use windows::write_pairwise_windows;

/// Format f64 for JSON: NaN/Infinity → null, -0 → 0, otherwise 6 decimal places.
fn format_json_f64(v: f64) -> String {
    if v.is_nan() || v.is_infinite() {
        "null".to_string()
    } else if v == 0.0 {
        "0.000000".to_string()
    } else {
        format!("{:.6}", v)
    }
}

/// Lineage plot data: (genome_id, lineage_name, dn_ds_ratio).
pub type LineagePlotResult = Vec<(String, String, f64)>;

/// Common output configuration shared across all write functions.
pub struct OutputConfig<'a> {
    /// Output file prefix (e.g., "results" → "results_pairwise_results.tsv")
    pub prefix: &'a str,
    /// Column separator character ('\\t' for TSV, ',' for CSV)
    pub sep: char,
    /// File extension ("tsv" or "csv")
    pub ext: &'a str,
    /// Which model was used (for column headers)
    pub model: Model,
    /// Optional summary statistics accumulator
    pub summary: Option<&'a SummaryStats>,
}

/// Spawn a writer thread that reassembles index-tagged blocks in ascending index
/// order, so the output is **deterministic** regardless of the order parallel tasks
/// finish (only the out-of-order window is buffered). Each task must send exactly one
/// `(index, block)` per index in `0..n`, including empty blocks so `index` advances.
///
/// The output file is created (and the header written) in the CALLING thread, so a
/// create/permission error surfaces as a clean `Err` instead of later panicking the
/// worker threads on a closed channel. The buffered writer is flushed explicitly so a
/// final-flush I/O error propagates rather than being swallowed by `Drop`.
type OrderedWriter = (
    crossbeam::channel::Sender<(usize, String)>,
    thread::JoinHandle<Result<(), std::io::Error>>,
);
fn spawn_ordered_writer(path: String, header: String) -> anyhow::Result<OrderedWriter> {
    let mut out = BufWriter::new(
        File::create(&path)
            .map_err(|e| anyhow::anyhow!("Cannot create '{}': {}", path, e))?,
    );
    out.write_all(header.as_bytes())?;
    let (tx, rx) = unbounded::<(usize, String)>();
    let handle = thread::spawn(move || -> Result<(), std::io::Error> {
        let mut pending: std::collections::HashMap<usize, String> = std::collections::HashMap::new();
        let mut next = 0usize;
        for (i, block) in rx {
            pending.insert(i, block);
            while let Some(b) = pending.remove(&next) {
                out.write_all(b.as_bytes())?;
                next += 1;
            }
        }
        while let Some(b) = pending.remove(&next) {
            out.write_all(b.as_bytes())?;
            next += 1;
        }
        out.flush()?;
        Ok(())
    });
    Ok((tx, handle))
}

