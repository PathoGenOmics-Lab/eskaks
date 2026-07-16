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

    // Load, validate, filter, and deduplicate sequences
    input::ensure_not_gzipped(&args.input_file)?;
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

    // Dispatch to the appropriate output mode
    let report_data = dispatch_output(&args, &data, &out_cfg, compute_pair, compute_pair_slices)?;

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
    }

    // Print summary
    if let Some(ref stats) = summary_stats {
        if args.summary {
            stats.print_summary();
        }
    }

    info!("All processes completed successfully.");
    Ok(())
}

/// Collect a fresh dN/dS distribution plus a capped `(dN, dS)` sample over
/// unique sequence pairs, for the report's always-present dN-vs-dS scatter and
/// histogram — independent of the primary output mode.
fn collect_report_pairwise(
    data: &input::SequenceData,
    engine: &ComputeEngine,
    scatter_cap: usize,
) -> (SummaryStats, Vec<(f64, f64)>) {
    let summary = SummaryStats::new();
    let n = data.n_unique;
    let total_pairs = n.saturating_sub(1) * n / 2;
    let stride = (total_pairs / scatter_cap.max(1)).max(1);
    let mut scatter = Vec::new();
    let mut local = stats::FloatAccum::new();
    let mut k = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            let (dn, ds) = engine.compute_slices(
                &data.unique_codon_indices[i],
                &data.unique_codon_indices[j],
            );
            let ratio = if ds > 0.0 {
                dn / ds
            } else if dn > 0.0 {
                f64::INFINITY
            } else {
                f64::NAN
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

/// Dispatch to the correct output mode based on CLI flags.
fn dispatch_output(
    args: &cli::FastaArgs,
    data: &input::SequenceData,
    cfg: &OutputConfig,
    compute_pair: impl Fn(usize, usize) -> DsDn + Sync,
    compute_pair_slices: impl Fn(&[u8], &[u8]) -> DsDn + Sync,
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
        info!("Results saved to {}_group_avg_dn_ds.{}", args.output, ext);

        if args.report {
            rd.group = Some(plot_data.clone());
        }
        if args.plot && !plot_data.is_empty() {
            let plot_path = format!("{}_group_dnds.svg", args.output);
            plot::group_bar_svg(&plot_data, &plot_path)?;
            info!("Plot saved to {}", plot_path);
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
        info!(
            "Lineage summary saved to {}_lineage_summary.{}",
            args.output, ext
        );

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
        }
    } else if let Some(win_size) = args.window_size {
        dispatch_window(args, data, cfg, compute_pair_slices, win_size)?;
    } else {
        info!("Generating pairwise results...");
        output::write_pairwise(
            &data.ids,
            &data.uidx_by_id,
            data.n_unique,
            compute_pair,
            cfg,
        )?;
        info!(
            "Results saved to {}_pairwise_results.{}",
            args.output, ext
        );

        if args.plot {
            if let Some(stats) = cfg.summary {
                let plot_path = format!("{}_dnds_histogram.svg", args.output);
                plot::histogram_svg(stats, &plot_path)?;
                info!("Plot saved to {}", plot_path);
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
    info!(
        "Results saved to {}_pairwise_windows.{}",
        args.output, cfg.ext
    );

    if let Some(ws) = &window_stats {
        let plot_path = format!("{}_window_plot.svg", args.output);
        plot::window_plot_svg(ws, &plot_path, win_size, args.window_step)?;
        info!("Plot saved to {}", plot_path);
    }
    Ok(())
}
