mod cli;
mod codon;
mod compute;
mod genetic_code;
mod gff;
mod input;
mod models;
mod output;
mod plot;
mod stats;
mod vcf;
mod vcf_analysis;

use anyhow::bail;
use clap::Parser;
use log::info;

use cli::{Args, SubCmd};
use compute::ComputeEngine;
use models::DsDn;
use output::OutputConfig;
use stats::SummaryStats;

fn main() -> anyhow::Result<()> {
    env_logger::init();
    let args = Args::parse();

    // Handle --list-codes (top-level flag)
    if args.list_codes {
        eprintln!("Available NCBI genetic code tables:");
        for (id, name) in genetic_code::list_tables() {
            eprintln!("  {:>2}  {}", id, name);
        }
        return Ok(());
    }

    match args.command {
        Some(SubCmd::Fasta(fasta_args)) => run_fasta(fasta_args),
        Some(SubCmd::Vcf(vcf_args)) => run_vcf(vcf_args),
        None => {
            // No subcommand: print help
            use clap::CommandFactory;
            Args::command().print_help()?;
            println!();
            Ok(())
        }
    }
}

fn run_fasta(args: cli::FastaArgs) -> anyhow::Result<()> {
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
    let stop_indices = genetic_code::stop_codon_indices(gc, args.model);
    let data = input::load_sequences(&args.input_file, args.model, args.min_codons, Some(&stop_indices))?;

    // Build compute engine
    let engine = ComputeEngine::new(args.model, gc);

    // Summary stats (only when --summary or --plot)
    let summary_stats = if args.summary || args.plot {
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
    dispatch_output(&args, &data, &out_cfg, compute_pair, compute_pair_slices)?;

    // Print summary
    if let Some(ref stats) = summary_stats {
        if args.summary {
            stats.print_summary();
        }
    }

    info!("All processes completed successfully.");
    Ok(())
}

fn run_vcf(args: cli::VcfArgs) -> anyhow::Result<()> {
    let gc = genetic_code::get_table(args.genetic_code).ok_or_else(|| {
        anyhow::anyhow!(
            "Unknown genetic code table {}. Use --list-codes to see available tables.",
            args.genetic_code
        )
    })?;
    if args.genetic_code != 1 {
        info!("Using genetic code table {}: {}", gc.id, gc.name);
    }

    let ref_path = std::path::Path::new(&args.reference);
    let gff_path = std::path::Path::new(&args.gff);
    let vcf_path = std::path::Path::new(&args.vcf);

    info!("Loading reference FASTA: {}", args.reference);
    let reference = vcf_analysis::parse_reference_fasta(ref_path)?;

    info!("Parsing GFF3 annotations: {}", args.gff);
    let genes = gff::parse_gff3(gff_path)?;
    info!("Found {} genes with CDS features", genes.len());

    info!("Parsing VCF file: {}", args.vcf);
    let snps = vcf::parse_vcf(vcf_path)?;
    info!("Found {} SNP records", snps.len());

    // Apply filters
    let snps = vcf::filter_snps(snps, args.pass_only, args.min_af, args.min_depth);
    info!("{} SNPs after filtering", snps.len());

    // Compute pN/pS
    let results = vcf_analysis::compute_pn_ps(&reference, &genes, &snps, gc);

    // Write results
    let output_path = vcf_analysis::write_results(&results, &args.output, &args.format)?;
    info!("Results saved to {}", output_path);

    // Generate plot if requested
    if args.plot {
        let plot_path = vcf_analysis::write_pnps_plot(&results, &args.output)?;
        info!("Plot saved to {}", plot_path);
    }

    // Print summary statistics
    let total_genes = results.len();
    let genes_with_snps = results.iter().filter(|r| r.total_snps > 0).count();
    let total_syn: u32 = results.iter().map(|r| r.syn_snps).sum();
    let total_nonsyn: u32 = results.iter().map(|r| r.nonsyn_snps).sum();
    eprintln!("\n── pN/pS Summary ──────────────────────────");
    eprintln!("  Genes analyzed:     {}", total_genes);
    eprintln!("  Genes with SNPs:    {}", genes_with_snps);
    eprintln!("  Total synonymous:   {}", total_syn);
    eprintln!("  Total nonsynonymous: {}", total_nonsyn);
    eprintln!("───────────────────────────────────────────");

    info!("VCF analysis completed successfully.");
    Ok(())
}

/// Dispatch to the correct output mode based on CLI flags.
fn dispatch_output(
    args: &cli::FastaArgs,
    data: &input::SequenceData,
    cfg: &OutputConfig,
    compute_pair: impl Fn(usize, usize) -> DsDn + Sync,
    compute_pair_slices: impl Fn(&[u8], &[u8]) -> DsDn + Sync,
) -> anyhow::Result<()> {
    let ext = cfg.ext;

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
    Ok(())
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
