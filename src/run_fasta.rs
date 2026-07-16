//! `eskaks fasta` subcommand: pairwise dN/dS orchestration and output dispatch.

use super::*;
use crate::compute::ComputeEngine;
use crate::output::OutputConfig;
use crate::models::DsDn;
use crate::stats::SummaryStats;
use anyhow::bail;
use log::info;

pub(crate) fn run_fasta(args: cli::FastaArgs) -> anyhow::Result<()> {
    // Validate genetic code
    let gc = genetic_code::get_table(args.genetic_code).ok_or_else(|| {
        anyhow::anyhow!(
            "Unknown genetic code table {}. Use --list-codes to see available tables.",
            args.genetic_code
        )
    })?;
    if args.genetic_code != 1 {
        info!("Using genetic code table {}: {}", gc.id, gc.name);
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.workers)
        .stack_size(4 * 1024 * 1024)
        .build_global()?;

    // Warn about flags that silently do nothing (or too little) on their own.
    if args.window_step != 1 && args.window_size.is_none() {
        log::warn!("--window-step is ignored without --window-size (no windowed analysis ran).");
    }
    if (1..100).contains(&args.bootstrap) {
        log::warn!(
            "--bootstrap {} is very low; a 95% CI from so few replicates is unreliable (use >= 1000).",
            args.bootstrap
        );
    }

    // Load, validate, filter, and deduplicate sequences. needletail decompresses
    // a gzipped FASTA transparently; we only guard against a directory here.
    input::ensure_not_directory(&args.input_file)?;
    let stop_indices = genetic_code::stop_codon_indices(gc, args.model);
    let data = input::load_sequences(&args.input_file, args.model, args.min_codons, Some(&stop_indices))?;

    // Build compute engine
    let engine = ComputeEngine::new(args.model, gc);

    // Summary stats (needed for --summary, --plot, and the --report histogram)
    let summary_stats = if args.summary || args.plot || args.report {
        Some(SummaryStats::new())
    } else {
        None
    };

    let ext = args.format.extension();
    let sep = args.format.separator();
    let out_cfg = OutputConfig {
        prefix: &args.output,
        sep,
        ext,
        model: args.model,
        summary: summary_stats.as_ref(),
    };

    // Closures that capture the engine
    let compute_pair = |u_i: usize, u_j: usize| -> DsDn { engine.compute_pair(&data, u_i, u_j) };

    let compute_pair_slices =
        |s1: &[u8], s2: &[u8]| -> DsDn {
            let (dn, ds) = engine.compute_slices(s1, s2);
            DsDn { dn, ds }
        };

    // Dispatch to the appropriate output mode, collecting every file written so the
    // run can confirm the output paths to the user (a plain run was otherwise silent).
    let mut written: Vec<String> = Vec::new();
    let report_data =
        dispatch_output(&args, &data, &out_cfg, compute_pair, compute_pair_slices, &mut written)?;

    // Interactive HTML report: a multi-panel dashboard. Always includes the
    // dN-vs-dS scatter and the dN/dS distribution (from a dedicated pairwise
    // pass), plus the lineage/group scatter when those modes were run.
    if args.report {
        let model_name = match args.model {
            models::Model::Nei => "Nei-Gojobori",
            models::Model::Li => "Li / LPB93",
        };
        let (rep_summary, dn_ds) = collect_report_pairwise(&data, &engine, 8000);
        let window = collect_window_profile(&data, &engine);
        let path = report::write_fasta_report(
            &args.output,
            model_name,
            Some(&rep_summary),
            report_data.lineage.as_deref(),
            report_data.group.as_deref(),
            Some(&dn_ds),
            Some(&window),
        )?;
        info!("Report saved to {}", path);
        written.push(path);
    }

    // Optional per-pair Nei-Gojobori neutrality test (variance + Z-test).
    if args.neutrality {
        if args.model != models::Model::Nei {
            info!("Note: analytic NG variances are Nei-only; SE/Z/P will be NaN for --model li.");
        }
        let stats = |u_i: usize, u_j: usize| engine.compute_pair_stats(&data, u_i, u_j);
        let path = output::write_pairwise_tests(
            &data.ids, &data.uidx_by_id, stats, &args.output, sep, ext,
        )?;
        info!("Neutrality test saved to {}", path);
        written.push(path);
    }

    // Optional per-pair bootstrap CIs (model-agnostic; resamples codons).
    if args.bootstrap > 0 {
        let seed = args.seed;
        let n_boot = args.bootstrap;
        let boot = |u_i: usize, u_j: usize| {
            let s1 = &data.unique_codon_indices[u_i];
            let s2 = &data.unique_codon_indices[u_j];
            let (dn, ds) = engine.compute_slices(s1, s2);
            // Per-pair seed keeps the bootstrap reproducible and independent
            // across pairs (so parallel writing stays deterministic).
            let pair_seed = seed
                ^ (u_i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
                ^ (u_j as u64).wrapping_mul(0xC2B2_AE3D_27D4_EB4F);
            let (dn_lo, dn_hi, ds_lo, ds_hi, r_lo, r_hi) =
                engine.bootstrap_ci(s1, s2, n_boot, pair_seed);
            (dn, ds, dn_lo, dn_hi, ds_lo, ds_hi, r_lo, r_hi)
        };
        let path = output::write_pairwise_bootstrap(
            &data.ids, &data.uidx_by_id, boot, &args.output, sep, ext,
        )?;
        info!("Bootstrap CIs saved to {}", path);
        written.push(path);
    }

    // Print summary
    if let Some(ref stats) = summary_stats {
        if args.summary {
            stats.print_summary();
        }
    }

    // Always confirm what was done and where the output went (unless silenced): a plain
    // `eskaks fasta in.fasta` previously finished with no terminal output at all.
    // `--quiet` drops the level to Error, so the Warn check doubles as the quiet gate.
    if log::log_enabled!(log::Level::Warn) {
        let model_name = match args.model {
            models::Model::Nei => "Nei-Gojobori",
            models::Model::Li => "Li / LPB93",
        };
        eprintln!("\n── Done ───────────────────────────────────");
        eprintln!(
            "  Sequences:  {} ({} unique) from {}",
            data.ids.len(), data.n_unique, args.input_file
        );
        eprintln!("  Model:      {}", model_name);
        if written.is_empty() {
            eprintln!("  Output:     (none)");
        } else {
            eprintln!("  Output:");
            for p in &written {
                eprintln!("    {}", p);
            }
        }
        eprintln!("────────────────────────────────────────────");
    }

    info!("All processes completed successfully.");
    Ok(())
}

/// Collect a fresh dN/dS distribution plus a capped `(dN, dS)` sample for the
/// report's always-present dN-vs-dS scatter and histogram — independent of the
/// primary output mode.
///
/// Iterates over ALL id-pairs (weighted by sequence multiplicity), not just unique
/// sequence pairs, so the report's Pairs / Valid pairs / Pooled / Mean cards match
/// the terminal `--summary` and the `_pairwise_results` table. Compute is still O(U)
/// per row via a per-row unique-index cache (mirroring `write_pairwise`), so duplicate
/// sequences — common in clonal datasets — are counted without recomputing.
fn collect_report_pairwise(
    data: &input::SequenceData,
    engine: &ComputeEngine,
    scatter_cap: usize,
) -> (SummaryStats, Vec<(f64, f64)>) {
    let summary = SummaryStats::new();
    let n_ids = data.ids.len();
    let n_u = data.n_unique;
    let total_pairs = n_ids.saturating_sub(1) * n_ids / 2;
    let stride = (total_pairs / scatter_cap.max(1)).max(1);
    let mut scatter = Vec::new();
    let mut local = stats::FloatAccum::new();
    let mut k = 0usize;
    let mut row_cache = vec![DsDn { dn: 0.0, ds: 0.0 }; n_u];
    let mut gen_map = vec![0u32; n_u];
    let mut cur_gen = 0u32;
    for i in 0..n_ids {
        cur_gen = cur_gen.wrapping_add(1);
        let u_i = data.uidx_by_id[i];
        row_cache[u_i] = engine.compute_pair(data, u_i, u_i);
        gen_map[u_i] = cur_gen;
        for j in (i + 1)..n_ids {
            let u_j = data.uidx_by_id[j];
            if gen_map[u_j] != cur_gen {
                row_cache[u_j] = engine.compute_pair(data, u_i, u_j);
                gen_map[u_j] = cur_gen;
            }
            let DsDn { dn, ds } = row_cache[u_j];
            // Use write_pairwise's exact ratio convention so the report histogram agrees
            // with --summary and the main table: a finite zero-divergence pair (dn == 0,
            // ds == 0 — e.g. identical/duplicate sequences) is 0.0, not NaN, so it bins
            // into [0.0, 0.2) rather than the [1.0, inf) overflow.
            let ratio = if dn.is_nan() || ds.is_nan() {
                f64::NAN
            } else if ds == 0.0 {
                if dn == 0.0 {
                    0.0
                } else {
                    f64::INFINITY
                }
            } else {
                dn / ds
            };
            summary.record_pair_atomic(dn, ds, ratio);
            local.record(dn, ds, ratio);
            if k.is_multiple_of(stride) && scatter.len() < scatter_cap {
                scatter.push((dn, ds));
            }
            k += 1;
        }
    }
    summary.flush_local(&local);
    (summary, scatter)
}

/// Sliding-window dN/dS profile along the alignment — a positional "Manhattan"
/// for the report. Returns `(codon_center, mean_dN/dS)` per window, or empty if
/// the sequences are not equal-length (unaligned) or too short. Uses a sampled
/// subset of pairs so it stays cheap.
fn collect_window_profile(data: &input::SequenceData, engine: &ComputeEngine) -> Vec<(usize, f64)> {
    let n = data.n_unique;
    if n < 2 {
        return Vec::new();
    }
    let seq_len = data.unique_codon_indices[0].len();
    if seq_len < 15 || data.unique_codon_indices.iter().any(|v| v.len() != seq_len) {
        return Vec::new();
    }
    let window = (seq_len / 15).clamp(10, seq_len);
    let step = (window / 2).max(1);

    // Sample up to ~300 pairs evenly across the pair space.
    let total = n.saturating_sub(1) * n / 2;
    let stride = (total / 300).max(1);
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    let mut k = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            if k.is_multiple_of(stride) {
                pairs.push((i, j));
            }
            k += 1;
        }
    }

    let mut profile = Vec::new();
    let mut start = 0;
    while start + window <= seq_len {
        let (mut sum, mut cnt) = (0.0f64, 0usize);
        for &(i, j) in &pairs {
            let (dn, ds) = engine.compute_slices(
                &data.unique_codon_indices[i][start..start + window],
                &data.unique_codon_indices[j][start..start + window],
            );
            if ds > 0.0 && dn.is_finite() {
                let r = dn / ds;
                if r.is_finite() {
                    sum += r;
                    cnt += 1;
                }
            }
        }
        if cnt > 0 {
            profile.push((start + window / 2, sum / cnt as f64));
        }
        start += step;
    }
    profile
}

/// Dispatch to the correct output mode based on CLI flags. Every file written is
/// pushed onto `written` so the caller can confirm the output paths to the user.
fn dispatch_output(
    args: &cli::FastaArgs,
    data: &input::SequenceData,
    cfg: &OutputConfig,
    compute_pair: impl Fn(usize, usize) -> DsDn + Sync,
    compute_pair_slices: impl Fn(&[u8], &[u8]) -> DsDn + Sync,
    written: &mut Vec<String>,
) -> anyhow::Result<report::FastaReportData> {
    let ext = cfg.ext;
    let mut rd = report::FastaReportData::default();

    if args.group_average {
        if args.window_size.is_some() {
            bail!("--window-size cannot be used with --group-average");
        }
        info!("Computing group average dN/dS...");
        let plot_data = output::write_group_average(
            &data.ids,
            &data.uidx_by_id,
            compute_pair,
            args.first_letter_lineage,
            cfg,
        )?;
        let out = format!("{}_group_avg_dn_ds.{}", args.output, ext);
        info!("Results saved to {}", out);
        written.push(out);

        if args.report {
            rd.group = Some(plot_data.clone());
        }
        if args.plot && !plot_data.is_empty() {
            let plot_path = format!("{}_group_dnds.svg", args.output);
            plot::group_bar_svg(&plot_data, &plot_path)?;
            info!("Plot saved to {}", plot_path);
            written.push(plot_path);
        }
    } else if args.lineage {
        if args.window_size.is_some() {
            bail!("--window-size cannot be used with --lineage");
        }
        info!("Computing dN/dS lineage summary...");
        let mut lineage_map: rustc_hash::FxHashMap<&str, usize> = rustc_hash::FxHashMap::default();
        let mut lineage_names: Vec<String> = Vec::new();
        let lineage_indices: Vec<usize> = data
            .ids
            .iter()
            .map(|id| {
                let key = codon::extract_group_key(id, args.first_letter_lineage);
                let next_idx = lineage_names.len();
                *lineage_map.entry(key).or_insert_with(|| {
                    lineage_names.push(key.to_string());
                    next_idx
                })
            })
            .collect();
        let plot_data = output::write_lineage(
            &data.ids,
            &data.uidx_by_id,
            data.n_unique,
            compute_pair,
            &lineage_indices,
            &lineage_names,
            cfg,
        )?;
        let out = format!("{}_lineage_summary.{}", args.output, ext);
        info!("Lineage summary saved to {}", out);
        written.push(out);

        if args.report {
            rd.lineage = Some(plot_data.clone());
        }
        if args.plot && !plot_data.is_empty() {
            let lineage_plot_data: Vec<plot::LineagePlotData> = plot_data
                .into_iter()
                .map(|(genome, lineage, ratio)| plot::LineagePlotData {
                    genome,
                    lineage,
                    ratio,
                })
                .collect();
            let plot_path = format!("{}_lineage_dnds.svg", args.output);
            plot::lineage_bar_svg(&lineage_plot_data, &plot_path)?;
            info!("Plot saved to {}", plot_path);
            written.push(plot_path);
        }
    } else if let Some(win_size) = args.window_size {
        dispatch_window(args, data, cfg, compute_pair_slices, win_size, written)?;
    } else {
        info!("Generating pairwise results...");
        output::write_pairwise(
            &data.ids,
            &data.uidx_by_id,
            data.n_unique,
            compute_pair,
            cfg,
        )?;
        let out = format!("{}_pairwise_results.{}", args.output, ext);
        info!("Results saved to {}", out);
        written.push(out);

        if args.plot {
            if let Some(stats) = cfg.summary {
                let plot_path = format!("{}_dnds_histogram.svg", args.output);
                plot::histogram_svg(stats, &plot_path)?;
                info!("Plot saved to {}", plot_path);
                written.push(plot_path);
            }
        }
    }
    Ok(rd)
}

/// Window mode dispatch with validation.
fn dispatch_window(
    args: &cli::FastaArgs,
    data: &input::SequenceData,
    cfg: &OutputConfig,
    compute_pair_slices: impl Fn(&[u8], &[u8]) -> DsDn + Sync,
    win_size: usize,
    written: &mut Vec<String>,
) -> anyhow::Result<()> {
    let seq_len = data
        .unique_codon_indices
        .first()
        .map(|v| v.len())
        .unwrap_or(0);
    if seq_len == 0 {
        bail!("Cannot use --window-size with empty sequences");
    }
    let misaligned = data
        .unique_codon_indices
        .iter()
        .any(|v| v.len() != seq_len);
    if misaligned {
        bail!("--window-size requires all sequences to have equal length. Sequences are not aligned");
    }
    if win_size == 0 || win_size > seq_len {
        bail!(
            "--window-size must be between 1 and {} (sequence length in codons)",
            seq_len
        );
    }
    if args.window_step == 0 {
        bail!("--window-step must be at least 1");
    }
    let num_windows = (seq_len - win_size) / args.window_step + 1;
    let window_stats = if args.plot {
        Some(stats::WindowStats::new(num_windows))
    } else {
        None
    };
    info!(
        "Generating sliding window pairwise results (window={}, step={})...",
        win_size, args.window_step
    );
    output::write_pairwise_windows(
        &data.ids,
        &data.uidx_by_id,
        &data.unique_codon_indices,
        compute_pair_slices,
        win_size,
        args.window_step,
        window_stats.as_ref(),
        cfg,
    )?;
    let out = format!("{}_pairwise_windows.{}", args.output, cfg.ext);
    info!("Results saved to {}", out);
    written.push(out);

    if let Some(ws) = &window_stats {
        let plot_path = format!("{}_window_plot.svg", args.output);
        plot::window_plot_svg(ws, &plot_path, win_size, args.window_step)?;
        info!("Plot saved to {}", plot_path);
        written.push(plot_path);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::SequenceData;
    use crate::models::Model;
    use std::sync::atomic::Ordering;

    #[test]
    fn report_pairwise_counts_all_id_pairs_not_unique() {
        // 4 ids but only 2 unique sequences (A, A, B, B). The report summary must count
        // all 4*3/2 = 6 id-pairs (weighted by sequence multiplicity), matching the
        // terminal --summary and the main table — not the single unique-sequence pair.
        let gc = crate::genetic_code::get_table(1).unwrap();
        let engine = ComputeEngine::new(Model::Nei, gc);
        let a = crate::codon::fasta_to_codon_indices(b"ATGGCTGCT", Model::Nei);
        let b = crate::codon::fasta_to_codon_indices(b"ATGATTGCT", Model::Nei);
        let data = SequenceData {
            ids: vec!["A0".into(), "A1".into(), "B0".into(), "B1".into()],
            uidx_by_id: vec![0, 0, 1, 1],
            unique_codon_indices: vec![a, b],
            n_unique: 2,
        };
        let (summary, _scatter) = collect_report_pairwise(&data, &engine, 8000);
        assert_eq!(
            summary.total_count.load(Ordering::Relaxed),
            6,
            "report must count all id-pairs (4*3/2 = 6), not unique-sequence pairs"
        );
    }

    #[test]
    fn report_histogram_bins_identical_pairs_like_the_main_table() {
        // Regression: identical (zero-divergence) pairs have dN/dS = 0.0 in the main
        // table and --summary (bin [0.0, 0.2)). The report's histogram must agree — it
        // must NOT map dn==0,ds==0 to NaN and land those pairs in the [1.0, inf) overflow
        // bin (which would read as positive selection on clonal data).
        let gc = crate::genetic_code::get_table(1).unwrap();
        let engine = ComputeEngine::new(Model::Nei, gc);
        let a = crate::codon::fasta_to_codon_indices(b"ATGGCTGCT", Model::Nei);
        let data = SequenceData {
            ids: vec!["s1".into(), "s2".into()], // identical -> one zero-divergence pair
            uidx_by_id: vec![0, 0],
            unique_codon_indices: vec![a],
            n_unique: 1,
        };
        let (summary, _scatter) = collect_report_pairwise(&data, &engine, 8000);
        let hist = summary.get_histogram();
        assert_eq!(hist[0].1, 1, "identical pair must bin into [0.0, 0.2), got {:?}", hist);
        assert_eq!(
            hist[hist.len() - 1].1,
            0,
            "identical pair must NOT land in the [1.0, inf) overflow bin, got {:?}",
            hist
        );
    }
}
