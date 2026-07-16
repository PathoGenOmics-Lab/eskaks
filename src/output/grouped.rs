//! Lineage and group-average table writers.

use super::*;

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

