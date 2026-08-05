//! Lineage and group-average table writers.

use super::*;

/// Format the `[low, high]` 95% CI field. In CSV mode the field is quoted so its
/// internal comma does not split it into two columns and misalign the row.
fn ci_field(lo: f64, hi: f64, sep: char) -> String {
    let inner = if lo.is_nan() && hi.is_nan() {
        "[NaN, NaN]".to_string()
    } else {
        format!("[{:.6}, {:.6}]", lo, hi)
    };
    if sep == ',' {
        format!("\"{}\"", inner)
    } else {
        inner
    }
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
    if num_lineages < 2 {
        log::warn!(
            "Only {} lineage detected from the sequence IDs. --lineage compares each genome ACROSS \
             lineages, so a single lineage gives a degenerate result — check the ID naming (grouping \
             splits on '_', or use --first-letter-lineage).",
            num_lineages
        );
    }
    let output_path = format!("{}_lineage_summary.{}", output_prefix, ext);
    let is_json = ext == "json";

    let (tx, writer_thread) = if is_json {
        spawn_ordered_json_writer(output_path.clone())?
    } else {
        spawn_ordered_writer(
            output_path.clone(),
            format!("Genome{s}Against_Lineage{s}Mean_dN{s}Mean_dS{s}dN/dS_Ratio\n", s = sep),
        )?
    };
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
                // Compute the self-comparison rather than hard-coding 0/0: identical
                // genomes deduplicate to this same unique index, so an all-N / all-gap
                // genome would otherwise contribute a spurious 0.0 to its lineage mean
                // instead of being excluded (compute_pair returns NaN for no comparable
                // codons). Mirrors the pairwise writer; compute_pair short-circuits u==u.
                row_cache[u_i] = compute_pair(u_i, u_i);
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
                    if is_json {
                        let _ = writeln!(block,
                            "{{\"genome\":\"{}\",\"against_lineage\":\"{}\",\"mean_dN\":{},\"mean_dS\":{},\"dN_dS\":{}}}",
                            json_escape(&ids[i]), json_escape(&lineage_names[lin_idx]),
                            format_json_f64(mean_dn), format_json_f64(mean_ds), format_json_f64(ratio));
                    } else {
                        let _ = writeln!(block, "{}{s}{}{s}{:.6}{s}{:.6}{s}{:.6}",
                            name_field(&ids[i], sep), name_field(&lineage_names[lin_idx], sep),
                            norm_zero(mean_dn), norm_zero(mean_ds), norm_zero(ratio), s = sep);
                    }

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
    if num_groups < 2 {
        log::warn!(
            "Only {} group detected from the sequence IDs. --group-average compares BETWEEN groups, \
             so a single group yields just a within-group row — check the ID naming (grouping splits \
             on '_', or use --first-letter-lineage).",
            num_groups
        );
    }
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
    let is_json = ext == "json";

    let (tx, writer_thread) = if is_json {
        spawn_ordered_json_writer(output_path.clone())?
    } else {
        spawn_ordered_writer(
            output_path.clone(),
            format!("Group1{s}Group2{s}NumSeqs1{s}NumSeqs2{s}NumComparisons{s}Mean_dN/dS{s}StdError{s}95%CI\n", s = sep),
        )?
    };
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
        // Local dN/dS accumulator so `--group-average --summary` reports the same dN/dS
        // min/max/mean and pooled lines that the pairwise and lineage summaries do.
        let mut local_stats = FloatAccum::new();

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
                    if let Some(stats) = summary {
                        stats.record_pair_atomic(result.dn, result.ds, ratio);
                        local_stats.record(result.dn, result.ds, ratio);
                    }
                } else if let Some(stats) = summary {
                    stats.record_pair_atomic(f64::NAN, f64::NAN, f64::NAN);
                }
            }
        }

        // Emit as a JSON object (NaN / N/A → null) or a delimited row. `mean`, `se`,
        // `ci_lo`, `ci_hi` are NaN where the corresponding column is N/A so the JSON and
        // plot paths share one source of truth.
        let g1n = &group_names[g1];
        let g2n = &group_names[g2];
        // Mean / SE / CI over the FINITE ratios only. A +inf ratio (dS = 0, dN > 0, i.e.
        // strong positive selection) is a real per-pair value but would poison an
        // arithmetic mean to inf and its variance to NaN, so it is excluded here exactly
        // as the pairwise summary excludes non-finite ratios. NumComparisons still counts
        // every defined pair; when no pair is finite the mean is N/A (NaN).
        let finite: Vec<f64> =
            pair_dn_ds_ratios.iter().copied().filter(|r| r.is_finite()).collect();
        let n_all = pair_dn_ds_ratios.len();
        let (n_comp, mean, se, ci_lo, ci_hi) = if finite.is_empty() {
            (n_all, f64::NAN, f64::NAN, f64::NAN, f64::NAN)
        } else {
            let n = finite.len();
            let mean: f64 = finite.iter().sum::<f64>() / n as f64;
            if n == 1 {
                (n_all, mean, f64::NAN, f64::NAN, f64::NAN)
            } else {
                let variance: f64 =
                    finite.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
                let se = (variance / n as f64).sqrt();
                let ci_hw = Z_95_CONFIDENCE * se;
                (n_all, mean, se, mean - ci_hw, mean + ci_hw)
            }
        };
        let plot_data = GroupPlotData {
            label: format!("{} vs {}", g1n, g2n),
            mean,
            ci_low: if n_comp == 1 { mean } else { ci_lo },
            ci_high: if n_comp == 1 { mean } else { ci_hi },
        };
        let line = if is_json {
            format!(
                "{{\"group1\":\"{}\",\"group2\":\"{}\",\"num_seqs1\":{},\"num_seqs2\":{},\"num_comparisons\":{},\"mean_dN_dS\":{},\"std_error\":{},\"ci_low\":{},\"ci_high\":{}}}\n",
                json_escape(g1n), json_escape(g2n), members1.len(), members2.len(), n_comp,
                format_json_f64(mean), format_json_f64(se), format_json_f64(ci_lo), format_json_f64(ci_hi)
            )
        } else if n_comp == 0 {
            format!("{}{s}{}{s}{}{s}{}{s}0{s}NaN{s}NaN{s}{}\n",
                name_field(g1n, sep), name_field(g2n, sep),
                members1.len(), members2.len(), ci_field(f64::NAN, f64::NAN, sep), s = sep)
        } else if n_comp == 1 {
            format!("{}{s}{}{s}{}{s}{}{s}{}{s}{:.6}{s}N/A{s}N/A\n",
                name_field(g1n, sep), name_field(g2n, sep),
                members1.len(), members2.len(), n_comp, mean, s = sep)
        } else {
            format!("{}{s}{}{s}{}{s}{}{s}{}{s}{:.6}{s}{:.6}{s}{}\n",
                name_field(g1n, sep), name_field(g2n, sep),
                members1.len(), members2.len(), n_comp,
                mean, se, ci_field(ci_lo, ci_hi, sep), s = sep)
        };

        if let Some(stats) = summary {
            stats.flush_local(&local_stats);
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

