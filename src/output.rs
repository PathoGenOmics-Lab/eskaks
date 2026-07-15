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

/// Writes pairwise results using a dedicated writer thread.
/// Computes pairs on-the-fly with lazy per-row caching (O(U) memory per thread).
/// Write a per-pair Nei-Gojobori neutrality-test table:
/// Seq1, Seq2, dN, dS, SE_dN, SE_dS, Z, P_value, where Z = (dN-dS)/√(V(dN)+V(dS))
/// and P_value is the two-sided normal p. `stats(u_i, u_j)` returns
/// `(dn, ds, var_dn, var_ds)`. SE/Z/P are NaN for the Li model (no NG variance).
pub fn write_pairwise_tests(
    ids: &[String],
    uidx_by_id: &[usize],
    stats: impl Fn(usize, usize) -> (f64, f64, f64, f64) + Sync,
    prefix: &str,
    sep: char,
    ext: &str,
) -> anyhow::Result<String> {
    use rayon::prelude::*;

    let n = ids.len();
    let pairs: Vec<(usize, usize)> =
        (0..n).flat_map(|i| (i + 1..n).map(move |j| (i, j))).collect();
    let is_json = ext == "json";
    let fmt = |v: f64| if v.is_finite() { format!("{:.6}", v) } else { "NaN".to_string() };
    let jfmt = |v: f64| if v.is_finite() { format!("{:.6}", v) } else { "null".to_string() };

    let rows: Vec<String> = pairs
        .par_iter()
        .map(|&(i, j)| {
            let (dn, ds, var_dn, var_ds) = stats(uidx_by_id[i], uidx_by_id[j]);
            let var_sum = var_dn + var_ds;
            let z = if var_sum > 0.0 { (dn - ds) / var_sum.sqrt() } else { f64::NAN };
            let p = crate::stats::normal_two_sided_p(z);
            if is_json {
                format!(
                    "  {{\"seq1\":\"{}\",\"seq2\":\"{}\",\"dN\":{},\"dS\":{},\"SE_dN\":{},\"SE_dS\":{},\"Z\":{},\"P_value\":{}}}",
                    ids[i], ids[j], jfmt(dn), jfmt(ds), jfmt(var_dn.sqrt()), jfmt(var_ds.sqrt()), jfmt(z), jfmt(p)
                )
            } else {
                format!(
                    "{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}",
                    ids[i], ids[j], fmt(dn), fmt(ds), fmt(var_dn.sqrt()), fmt(var_ds.sqrt()), fmt(z), fmt(p), s = sep
                )
            }
        })
        .collect();

    let output_path = format!("{}_pairwise_tests.{}", prefix, ext);
    let mut file = BufWriter::new(File::create(&output_path)?);
    if is_json {
        writeln!(file, "[")?;
        writeln!(file, "{}", rows.join(",\n"))?;
        writeln!(file, "]")?;
    } else {
        writeln!(file, "Seq1{s}Seq2{s}dN{s}dS{s}SE_dN{s}SE_dS{s}Z{s}P_value", s = sep)?;
        for row in &rows {
            writeln!(file, "{}", row)?;
        }
    }
    Ok(output_path)
}

/// Write a per-pair bootstrap-CI table:
/// Seq1, Seq2, dN, dN_CI_low, dN_CI_high, dS, dS_CI_low, dS_CI_high,
/// dN/dS, dNdS_CI_low, dNdS_CI_high. `stats(u_i, u_j)` returns
/// `(dn, ds, dn_lo, dn_hi, ds_lo, ds_hi, ratio_lo, ratio_hi)`.
pub fn write_pairwise_bootstrap(
    ids: &[String],
    uidx_by_id: &[usize],
    stats: impl Fn(usize, usize) -> (f64, f64, f64, f64, f64, f64, f64, f64) + Sync,
    prefix: &str,
    sep: char,
    ext: &str,
) -> anyhow::Result<String> {
    use rayon::prelude::*;

    let n = ids.len();
    let pairs: Vec<(usize, usize)> =
        (0..n).flat_map(|i| (i + 1..n).map(move |j| (i, j))).collect();
    let is_json = ext == "json";
    let fmt = |v: f64| if v.is_finite() { format!("{:.6}", v) } else { "NaN".to_string() };
    let jfmt = |v: f64| if v.is_finite() { format!("{:.6}", v) } else { "null".to_string() };

    let rows: Vec<String> = pairs
        .par_iter()
        .map(|&(i, j)| {
            let (dn, ds, dn_lo, dn_hi, ds_lo, ds_hi, r_lo, r_hi) = stats(uidx_by_id[i], uidx_by_id[j]);
            let ratio = if ds > 0.0 { dn / ds } else { f64::NAN };
            if is_json {
                format!(
                    "  {{\"seq1\":\"{}\",\"seq2\":\"{}\",\"dN\":{},\"dN_CI_low\":{},\"dN_CI_high\":{},\"dS\":{},\"dS_CI_low\":{},\"dS_CI_high\":{},\"dN_dS\":{},\"dNdS_CI_low\":{},\"dNdS_CI_high\":{}}}",
                    ids[i], ids[j], jfmt(dn), jfmt(dn_lo), jfmt(dn_hi), jfmt(ds), jfmt(ds_lo), jfmt(ds_hi), jfmt(ratio), jfmt(r_lo), jfmt(r_hi)
                )
            } else {
                format!(
                    "{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}",
                    ids[i], ids[j], fmt(dn), fmt(dn_lo), fmt(dn_hi), fmt(ds), fmt(ds_lo), fmt(ds_hi), fmt(ratio), fmt(r_lo), fmt(r_hi), s = sep
                )
            }
        })
        .collect();

    let output_path = format!("{}_pairwise_bootstrap.{}", prefix, ext);
    let mut file = BufWriter::new(File::create(&output_path)?);
    if is_json {
        writeln!(file, "[")?;
        writeln!(file, "{}", rows.join(",\n"))?;
        writeln!(file, "]")?;
    } else {
        writeln!(file, "Seq1{s}Seq2{s}dN{s}dN_CI_low{s}dN_CI_high{s}dS{s}dS_CI_low{s}dS_CI_high{s}dN/dS{s}dNdS_CI_low{s}dNdS_CI_high", s = sep)?;
        for row in &rows {
            writeln!(file, "{}", row)?;
        }
    }
    Ok(output_path)
}

pub fn write_pairwise(
    ids: &[String],
    uidx_by_id: &[usize],
    n_u: usize,
    compute_pair: impl Fn(usize, usize) -> DsDn + Sync,
    cfg: &OutputConfig,
) -> anyhow::Result<()> {
    let output_prefix = cfg.prefix;
    let model = cfg.model;
    let sep = cfg.sep;
    let ext = cfg.ext;
    let summary = cfg.summary;
    let total_pairs_to_write = ids.len() * (ids.len() - 1) / 2;
    let pb_write = ProgressBar::new(total_pairs_to_write as u64);
    pb_write.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] Computing & writing pairs: {pos}/{len} ({eta})")
        .progress_chars("#>-"));

    let nan_count = AtomicUsize::new(0);
    // Each parallel task sends exactly one (row_index, block) tuple. The writer
    // reassembles blocks in ascending row order so the output is deterministic
    // regardless of thread scheduling (only the out-of-order window is buffered).
    let (tx, rx) = unbounded::<(usize, String)>();

    let is_json = ext == "json";
    // Create the output file in this (parent) thread so a create/permission error is
    // returned cleanly instead of panicking the workers on a closed channel later.
    let out_path = format!("{}_pairwise_results.{}", output_prefix, ext);
    let mut out_file = BufWriter::new(
        File::create(&out_path).map_err(|e| anyhow::anyhow!("Cannot create '{}': {}", out_path, e))?,
    );
    let writer_thread = thread::spawn({
        move || -> Result<(), std::io::Error> {
            let mut pending: std::collections::HashMap<usize, String> = std::collections::HashMap::new();
            let mut next = 0usize;
            if is_json {
                out_file.write_all(b"[\n")?;
                let mut first = true;
                let emit = |out_file: &mut BufWriter<File>, first: &mut bool, block: &str| -> Result<(), std::io::Error> {
                    for line in block.lines() {
                        if !*first { out_file.write_all(b",\n")?; }
                        out_file.write_all(line.as_bytes())?;
                        *first = false;
                    }
                    Ok(())
                };
                for (i, block) in rx {
                    pending.insert(i, block);
                    while let Some(b) = pending.remove(&next) { emit(&mut out_file, &mut first, &b)?; next += 1; }
                }
                while let Some(b) = pending.remove(&next) { emit(&mut out_file, &mut first, &b)?; next += 1; }
                out_file.write_all(b"\n]\n")?;
            } else {
                let header = match model {
                    Model::Li => format!("Seq1{s}Seq2{s}dN(Ka){s}dS(Ks){s}dN/dS\n", s = sep),
                    Model::Nei => format!("Seq1{s}Seq2{s}dN{s}dS{s}dN/dS\n", s = sep),
                };
                out_file.write_all(header.as_bytes())?;
                for (i, block) in rx {
                    pending.insert(i, block);
                    while let Some(b) = pending.remove(&next) { out_file.write_all(b.as_bytes())?; next += 1; }
                }
                while let Some(b) = pending.remove(&next) { out_file.write_all(b.as_bytes())?; next += 1; }
            }
            out_file.flush()?;
            Ok(())
        }
    });

    (0..ids.len()).into_par_iter().for_each_init(
        || {
            (
                tx.clone(),
                vec![DsDn { dn: 0.0, ds: 0.0 }; n_u],
                vec![0u32; n_u],
                0u32,
                String::with_capacity(1024 * 64),
                FloatAccum::new(),
            )
        },
        |(sender, row_cache, gen_map, cur_gen, local_buffer, local_stats), i| {
            *cur_gen = cur_gen.wrapping_add(1);
            let gen = *cur_gen;
            let u_i = uidx_by_id[i];
            row_cache[u_i] = DsDn { dn: 0.0, ds: 0.0 };
            gen_map[u_i] = gen;

            for j in (i + 1)..ids.len() {
                let u_j = uidx_by_id[j];
                if gen_map[u_j] != gen {
                    row_cache[u_j] = compute_pair(u_i, u_j);
                    gen_map[u_j] = gen;
                }
                let result = row_cache[u_j];
                if !result.dn.is_finite() || !result.ds.is_finite() {
                    nan_count.fetch_add(1, Ordering::Relaxed);
                }
                let ratio = if result.ds == 0.0 {
                    if result.dn == 0.0 { 0.0 } else { f64::INFINITY }
                } else {
                    result.dn / result.ds
                };

                if let Some(stats) = summary {
                    stats.record_pair_atomic(result.dn, result.ds, ratio);
                    local_stats.record(result.dn, result.ds, ratio);
                }

                if is_json {
                    let _ = writeln!(local_buffer,
                        "{{\"seq1\":\"{}\",\"seq2\":\"{}\",\"dN\":{},\"dS\":{},\"dN_dS\":{}}}",
                        &ids[i], &ids[j],
                        format_json_f64(result.dn), format_json_f64(result.ds), format_json_f64(ratio));
                } else {
                    let _ = writeln!(local_buffer, "{}{s}{}{s}{:.6}{s}{:.6}{s}{:.6}",
                        &ids[i], &ids[j], result.dn, result.ds, ratio, s = sep);
                }

            }

            // Flush thread-local stats to shared accumulator (once per row)
            if let Some(stats) = summary {
                stats.flush_local(local_stats);
                local_stats.reset();
            }

            // Send this row's block exactly once, tagged with its index, so the
            // writer can emit rows in order (empty for the last row, which has no pairs).
            sender.send((i, std::mem::take(local_buffer)))
                .expect("Writer thread channel closed unexpectedly");
            pb_write.inc((ids.len() - 1 - i) as u64);
        },
    );
    drop(tx);

    writer_thread.join()
        .expect("Writer thread panicked")
        .expect("Writer thread encountered an I/O error");
    pb_write.finish_with_message("Pairwise computation & writing completed.");

    let nan_total = nan_count.load(Ordering::Relaxed);
    if nan_total > 0 {
        info!("{} of {} pairs ({:.1}%) returned NaN due to saturation.",
            nan_total, total_pairs_to_write,
            100.0 * nan_total as f64 / total_pairs_to_write as f64);
    }
    Ok(())
}

/// Writes dN/dS summary by lineage using a dedicated writer thread.
/// Computes pairs on-the-fly with lazy per-row caching.
/// Returns lineage plot data if summary stats are being collected.
pub fn write_lineage(
    ids: &[String],
    uidx_by_id: &[usize],
    n_u: usize,
    compute_pair: impl Fn(usize, usize) -> DsDn + Sync,
    lineage_indices: &[usize],
    lineage_names: &[String],
    cfg: &OutputConfig,
) -> anyhow::Result<LineagePlotResult> {
    let output_prefix = cfg.prefix;
    let sep = cfg.sep;
    let ext = cfg.ext;
    let summary = cfg.summary;
    let num_lineages = lineage_names.len();
    let output_path = format!("{}_lineage_summary.{}", output_prefix, ext);

    let (tx, writer_thread) = spawn_ordered_writer(
        output_path.clone(),
        format!("Genome{s}Against_Lineage{s}Mean_dN{s}Mean_dS{s}dN/dS_Ratio\n", s = sep),
    )?;
    let (plot_tx, plot_rx) = if summary.is_some() {
        let (t, r) = unbounded::<(String, String, f64)>();
        (Some(t), Some(r))
    } else {
        (None, None)
    };

    (0..ids.len())
        .into_par_iter()
        .for_each_init(
            || {
                (
                    tx.clone(),
                    plot_tx.clone(),
                    vec![DsDn { dn: 0.0, ds: 0.0 }; n_u],
                    vec![0u32; n_u],
                    0u32,
                    vec![(0.0f64, 0.0f64, 0usize); num_lineages],
                    FloatAccum::new(),
                )
            },
            |(sender, plot_sender, row_cache, gen_map, cur_gen, local_aggr, local_stats), i| {
                *cur_gen = cur_gen.wrapping_add(1);
                let gen = *cur_gen;
                let u_i = uidx_by_id[i];
                row_cache[u_i] = DsDn { dn: 0.0, ds: 0.0 };
                gen_map[u_i] = gen;

                for a in local_aggr.iter_mut() { *a = (0.0, 0.0, 0); }

                for j in 0..ids.len() {
                    if i == j { continue; }
                    let u_j = uidx_by_id[j];
                    if gen_map[u_j] != gen {
                        row_cache[u_j] = compute_pair(u_i, u_j);
                        gen_map[u_j] = gen;
                    }
                    let result = row_cache[u_j];

                    if result.dn.is_finite() && result.ds.is_finite() {
                        let lin_idx = lineage_indices[j];
                        local_aggr[lin_idx].0 += result.dn;
                        local_aggr[lin_idx].1 += result.ds;
                        local_aggr[lin_idx].2 += 1;
                    }
                }

                let mut block = String::new();
                for (lin_idx, &(sum_dn, sum_ds, count)) in local_aggr.iter().enumerate() {
                    if count == 0 { continue; }
                    let mean_dn = sum_dn / count as f64;
                    let mean_ds = sum_ds / count as f64;
                    let ratio = if mean_ds == 0.0 {
                        if mean_dn == 0.0 { 0.0 } else { f64::INFINITY }
                    } else {
                        mean_dn / mean_ds
                    };
                    let _ = writeln!(block, "{}{s}{}{s}{:.6}{s}{:.6}{s}{:.6}",
                        &ids[i], &lineage_names[lin_idx], mean_dn, mean_ds, ratio, s = sep);

                    if let Some(stats) = summary {
                        stats.record_pair_atomic(mean_dn, mean_ds, ratio);
                        local_stats.record(mean_dn, mean_ds, ratio);
                    }
                    if let Some(ref ps) = plot_sender {
                        let _ = ps.send((ids[i].clone(), lineage_names[lin_idx].clone(), ratio));
                    }
                }

                if let Some(stats) = summary {
                    stats.flush_local(local_stats);
                    local_stats.reset();
                }

                // Send once per row (even if empty) so the ordered writer can advance.
                sender.send((i, block)).expect("Writer thread channel closed unexpectedly");
            },
        );
    drop(tx);
    drop(plot_tx);

    writer_thread.join()
        .expect("Writer thread panicked")
        .expect("Writer thread encountered an I/O error");

    let plot_data = if let Some(rx) = plot_rx {
        // Received in thread-scheduling order; sort for a deterministic plot.
        let mut v: Vec<(String, String, f64)> = rx.into_iter().collect();
        v.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
        v
    } else {
        Vec::new()
    };

    Ok(plot_data)
}

/// Writes grouped dN/dS averages using a dedicated writer thread.
/// Returns group plot data if summary stats are being collected.
pub fn write_group_average(
    ids: &[String],
    uidx_by_id: &[usize],
    compute_pair: impl Fn(usize, usize) -> DsDn + Sync,
    first_letter_lineage: bool,
    cfg: &OutputConfig,
) -> anyhow::Result<Vec<GroupPlotData>> {
    let output_prefix = cfg.prefix;
    let sep = cfg.sep;
    let ext = cfg.ext;
    let summary = cfg.summary;
    let mut group_map: rustc_hash::FxHashMap<&str, usize> = rustc_hash::FxHashMap::default();
    let mut group_names: Vec<String> = Vec::new();
    let group_by_id: Vec<usize> = ids.iter().map(|id| {
        let key = extract_group_key(id, first_letter_lineage);
        let next_idx = group_names.len();
        *group_map.entry(key).or_insert_with(|| {
            group_names.push(key.to_string());
            next_idx
        })
    }).collect();

    let num_groups = group_names.len();
    let mut group_members: Vec<Vec<usize>> = vec![Vec::new(); num_groups];
    for (id_idx, &grp) in group_by_id.iter().enumerate() {
        group_members[grp].push(id_idx);
    }

    let mut group_pairs: Vec<(usize, usize)> = Vec::new();
    for i in 0..num_groups {
        for j in i..num_groups {
            group_pairs.push((i, j));
        }
    }

    let pb_group = ProgressBar::new(group_pairs.len() as u64);
    pb_group.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] Groups: {pos}/{len} ({eta})")
            .progress_chars("#>-"),
    );

    let output_path = format!("{}_group_avg_dn_ds.{}", output_prefix, ext);

    let (tx, writer_thread) = spawn_ordered_writer(
        output_path.clone(),
        format!("Group1{s}Group2{s}NumSeqs1{s}NumSeqs2{s}NumComparisons{s}Mean_dN/dS{s}StdError{s}95%CI\n", s = sep),
    )?;
    let (plot_tx, plot_rx) = if summary.is_some() {
        let (t, r) = unbounded::<GroupPlotData>();
        (Some(t), Some(r))
    } else {
        (None, None)
    };

    group_pairs.into_par_iter().enumerate().for_each_with((tx, plot_tx), |(s, ps), (pair_idx, (g1, g2))| {
        let members1 = &group_members[g1];
        let members2 = &group_members[g2];
        let mut pair_dn_ds_ratios = Vec::new();
        let mut nan_pair_count = 0usize;

        for (pos, &id1_idx) in members1.iter().enumerate() {
            let u1 = uidx_by_id[id1_idx];
            let start = if g1 == g2 { pos + 1 } else { 0 };
            for &id2_idx in &members2[start..] {
                let u2 = uidx_by_id[id2_idx];
                let result = compute_pair(u1, u2);
                if result.dn.is_finite() && result.ds.is_finite() {
                    let ratio = if result.ds == 0.0 {
                        if result.dn == 0.0 { 0.0 } else { f64::INFINITY }
                    } else {
                        result.dn / result.ds
                    };
                    pair_dn_ds_ratios.push(ratio);
                } else {
                    nan_pair_count += 1;
                }
            }
        }

        let (line, plot_data) = if pair_dn_ds_ratios.is_empty() {
            (format!("{}{s}{}{s}{}{s}{}{s}0{s}NaN{s}NaN{s}[NaN, NaN]\n",
                &group_names[g1], &group_names[g2],
                members1.len(), members2.len(), s = sep),
             GroupPlotData {
                 label: format!("{} vs {}", &group_names[g1], &group_names[g2]),
                 mean: f64::NAN, ci_low: f64::NAN, ci_high: f64::NAN,
             })
        } else {
            let n = pair_dn_ds_ratios.len();
            let mean: f64 = pair_dn_ds_ratios.iter().sum::<f64>() / n as f64;
            if n == 1 {
                (format!("{}{s}{}{s}{}{s}{}{s}{}{s}{:.6}{s}N/A{s}N/A\n",
                    &group_names[g1], &group_names[g2],
                    members1.len(), members2.len(), n, mean, s = sep),
                 GroupPlotData {
                     label: format!("{} vs {}", &group_names[g1], &group_names[g2]),
                     mean, ci_low: mean, ci_high: mean,
                 })
            } else {
                let variance: f64 =
                    pair_dn_ds_ratios.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
                let se = (variance / n as f64).sqrt();
                let ci_hw = Z_95_CONFIDENCE * se;
                (format!("{}{s}{}{s}{}{s}{}{s}{}{s}{:.6}{s}{:.6}{s}[{:.6}, {:.6}]\n",
                    &group_names[g1], &group_names[g2],
                    members1.len(), members2.len(), n,
                    mean, se, mean - ci_hw, mean + ci_hw, s = sep),
                 GroupPlotData {
                     label: format!("{} vs {}", &group_names[g1], &group_names[g2]),
                     mean, ci_low: mean - ci_hw, ci_high: mean + ci_hw,
                 })
            }
        };

        if let Some(stats) = summary {
            for _ in 0..nan_pair_count {
                stats.record_pair_atomic(f64::NAN, f64::NAN, f64::NAN);
            }
            for &r in &pair_dn_ds_ratios {
                stats.record_pair_atomic(0.0, 0.0, r);
            }
        }

        if let Some(ref plot_s) = ps {
            let _ = plot_s.send(plot_data);
        }

        s.send((pair_idx, line)).expect("Writer thread channel closed unexpectedly");
        pb_group.inc(1);
    });

    writer_thread.join()
        .expect("Writer thread panicked")
        .expect("Writer thread encountered an I/O error");
    pb_group.finish_with_message("Group average computation completed.");

    let plot_data = if let Some(rx) = plot_rx {
        // Received in thread-scheduling order; sort by label for a deterministic plot.
        let mut v: Vec<GroupPlotData> = rx.into_iter().collect();
        v.sort_by(|a, b| a.label.cmp(&b.label));
        v
    } else {
        Vec::new()
    };

    Ok(plot_data)
}

/// Writes pairwise sliding window results using a dedicated writer thread.
#[allow(clippy::too_many_arguments)]
pub fn write_pairwise_windows(
    ids: &[String],
    uidx_by_id: &[usize],
    unique_codon_indices: &[Vec<u8>],
    compute_pair_slices: impl Fn(&[u8], &[u8]) -> DsDn + Sync,
    window_size: usize,
    window_step: usize,
    window_stats: Option<&WindowStats>,
    cfg: &OutputConfig,
) -> anyhow::Result<()> {
    let output_prefix = cfg.prefix;
    let model = cfg.model;
    let sep = cfg.sep;
    let ext = cfg.ext;
    let summary = cfg.summary;
    let seq_len = unique_codon_indices.first().map(|v| v.len()).unwrap_or(0);
    let num_windows = if seq_len >= window_size {
        (seq_len - window_size) / window_step + 1
    } else {
        0
    };
    let total_pairs = ids.len() * (ids.len() - 1) / 2;
    let pb = ProgressBar::new((total_pairs * num_windows) as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] Window pairs: {pos}/{len} ({eta})")
        .progress_chars("#>-"));

    let nan_count = AtomicUsize::new(0);
    let header = match model {
        Model::Li => format!("Seq1{s}Seq2{s}Window_Start{s}Window_End{s}dN(Ka){s}dS(Ks){s}dN/dS\n", s = sep),
        Model::Nei => format!("Seq1{s}Seq2{s}Window_Start{s}Window_End{s}dN{s}dS{s}dN/dS\n", s = sep),
    };
    let (tx, writer_thread) =
        spawn_ordered_writer(format!("{}_pairwise_windows.{}", output_prefix, ext), header)?;

    (0..ids.len()).into_par_iter().for_each_init(
        || {
            (
                tx.clone(),
                String::with_capacity(1024 * 64),
                FloatAccum::new(),
            )
        },
        |(sender, local_buffer, local_stats), i| {
            let u_i = uidx_by_id[i];
            for j in (i + 1)..ids.len() {
                let u_j = uidx_by_id[j];
                let s1 = &unique_codon_indices[u_i];
                let s2 = &unique_codon_indices[u_j];
                for w in 0..num_windows {
                    let start = w * window_step;
                    let end = start + window_size;
                    let result = if u_i == u_j {
                        DsDn { dn: 0.0, ds: 0.0 }
                    } else {
                        compute_pair_slices(&s1[start..end], &s2[start..end])
                    };
                    if !result.dn.is_finite() || !result.ds.is_finite() {
                        nan_count.fetch_add(1, Ordering::Relaxed);
                    }
                    let ratio = if result.ds == 0.0 {
                        if result.dn == 0.0 { 0.0 } else { f64::INFINITY }
                    } else {
                        result.dn / result.ds
                    };

                    if let Some(stats) = summary {
                        stats.record_pair_atomic(result.dn, result.ds, ratio);
                        local_stats.record(result.dn, result.ds, ratio);
                    }
                    if let Some(ws) = window_stats {
                        ws.record(w, ratio);
                    }

                    let _ = writeln!(local_buffer, "{}{s}{}{s}{}{s}{}{s}{:.6}{s}{:.6}{s}{:.6}",
                        &ids[i], &ids[j], start + 1, end, result.dn, result.ds, ratio, s = sep);
                }
                pb.inc(num_windows as u64);
            }

            if let Some(stats) = summary {
                stats.flush_local(local_stats);
                local_stats.reset();
            }

            // One block per row i, in order, so the output is deterministic.
            sender.send((i, std::mem::take(local_buffer)))
                .expect("Writer thread channel closed unexpectedly");
        },
    );
    drop(tx);

    writer_thread.join()
        .expect("Writer thread panicked")
        .expect("Writer thread encountered an I/O error");
    pb.finish_with_message("Sliding window computation completed.");

    let nan_total = nan_count.load(Ordering::Relaxed);
    let total_computations = total_pairs * num_windows;
    if nan_total > 0 {
        info!("{} of {} window pairs ({:.1}%) returned NaN due to saturation.",
            nan_total, total_computations,
            100.0 * nan_total as f64 / total_computations as f64);
    }
    Ok(())
}
