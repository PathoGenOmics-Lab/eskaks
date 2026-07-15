mod cli;
mod codon;
mod compute;
mod genetic_code;
mod gff;
mod input;
mod models;
mod output;
mod plot;
mod report;
mod stats;
mod vcf;
mod vcf_analysis;

use anyhow::{bail, Context};
use clap::Parser;
use log::{info, warn, LevelFilter};

use std::collections::HashSet;

use cli::{Args, SubCmd};
use compute::ComputeEngine;
use models::DsDn;
use output::OutputConfig;
use stats::SummaryStats;

/// Format up to three names from a set for a diagnostic message.
fn sample_names(set: &HashSet<&str>) -> String {
    let mut names: Vec<&str> = set.iter().copied().collect();
    names.sort_unstable();
    let shown = names.len().min(3);
    let mut out = names[..shown].join(", ");
    if names.len() > shown {
        out.push_str(&format!(", … (+{} more)", names.len() - shown));
    }
    out
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    // Show data-quality warnings by default (previously env_logger defaulted to
    // "off", hiding every REF-mismatch / skipped-gene / saturation diagnostic
    // unless the user knew to set RUST_LOG). RUST_LOG still overrides this.
    let level = if args.quiet {
        LevelFilter::Error
    } else {
        match args.verbose {
            0 => LevelFilter::Warn,
            1 => LevelFilter::Info,
            _ => LevelFilter::Debug,
        }
    };
    env_logger::Builder::new()
        .filter_level(level)
        .parse_default_env()
        .init();

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

    if !args.kappa.is_finite() || args.kappa <= 0.0 {
        bail!("--kappa must be a positive, finite number (got {})", args.kappa);
    }
    if args.kappa != 1.0 {
        info!("Using mutation-spectrum-weighted site counting (kappa = {})", args.kappa);
    }

    // Validate allele-frequency filter ranges up front (a fat-fingered
    // --min-af 30 instead of 0.30 otherwise silently filters out every SNP).
    for (name, v) in [("--min-af", args.min_af), ("--max-af", args.max_af)] {
        if let Some(af) = v {
            if !(0.0..=1.0).contains(&af) {
                bail!("{} must be between 0.0 and 1.0 (got {})", name, af);
            }
        }
    }
    if let (Some(min), Some(max)) = (args.min_af, args.max_af) {
        if min > max {
            bail!("--min-af ({}) must not exceed --max-af ({})", min, max);
        }
    }
    if !args.fdr.is_finite() || args.fdr <= 0.0 || args.fdr > 1.0 {
        bail!("--fdr must be in (0, 1] (got {})", args.fdr);
    }
    if !(0.0..=1.0).contains(&args.mk_fixed_af) {
        bail!("--mk-fixed-af must be between 0.0 and 1.0 (got {})", args.mk_fixed_af);
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.workers)
        .stack_size(4 * 1024 * 1024)
        .build_global()?;

    let ref_path = std::path::Path::new(&args.reference);
    let gff_path = std::path::Path::new(&args.gff);

    info!("Loading reference FASTA: {}", args.reference);
    let reference = vcf_analysis::parse_reference_fasta(ref_path)?;

    info!("Parsing GFF3 annotations: {}", args.gff);
    let genes = gff::parse_gff3(gff_path)?;
    info!("Found {} genes with CDS features", genes.len());

    // Fail early on a total contig-name mismatch between GFF and reference —
    // otherwise every gene is silently skipped and the run "succeeds" with an
    // all-NaN output file (the classic 'Chromosome' vs 'NC_000962.3' footgun).
    {
        let ref_ids: HashSet<&str> = reference.keys().map(String::as_str).collect();
        let gene_ids: HashSet<&str> = genes.iter().map(|g| g.seqid.as_str()).collect();
        if !gene_ids.is_empty() && gene_ids.is_disjoint(&ref_ids) {
            bail!(
                "No GFF sequence name matches the reference FASTA.\n  GFF uses e.g.:       {}\n  reference has:       {}\nContig names must match across reference, GFF, and VCF.",
                sample_names(&gene_ids),
                sample_names(&ref_ids)
            );
        }
    }

    // Collect all VCF paths
    let mut vcf_paths: Vec<String> = args.vcf.clone();
    if let Some(ref list_path) = args.vcf_list {
        let content = std::fs::read_to_string(list_path)
            .with_context(|| format!("Failed to read VCF list file: {}", list_path))?;
        for line in content.lines() {
            let line = line.trim();
            if !line.is_empty() && !line.starts_with('#') {
                vcf_paths.push(line.to_string());
            }
        }
    }
    if vcf_paths.is_empty() {
        bail!("No VCF files provided. Use --vcf <file> or --vcf-list <file>");
    }

    // Parse and merge SNPs from all VCF files.
    // When multiple VCFs (one per sample), AF = fraction of samples carrying the ALT.
    let n_samples = vcf_paths.len();
    info!("Loading {} VCF file(s)...", n_samples);

    let snps = if n_samples == 1 {
        // Single VCF: use AF as-is (could be multi-sample VCF)
        let snps = vcf::parse_vcf(std::path::Path::new(&vcf_paths[0]))?;
        info!("Found {} SNP records", snps.len());
        vcf::filter_snps(snps, args.pass_only, args.min_af, args.max_af, args.min_depth)
    } else {
        // Multiple single-sample VCFs: merge and compute AF as fraction of samples
        let merged = vcf::merge_vcfs(&vcf_paths, args.pass_only, args.min_depth)?;
        info!("Merged {} unique SNP positions from {} samples", merged.len(), n_samples);
        vcf::filter_snps(merged, false, args.min_af, args.max_af, None)
    };
    info!("{} SNPs after filtering", snps.len());

    if snps.is_empty() {
        // Not fatal (0 SNPs is a valid, if uninformative, result), but make it
        // loud so it is never mistaken for a silently-empty output file.
        warn!(
            "0 SNPs remain after parsing and filtering — every gene will report \
             0 SNPs and undefined pN/pS. Check the VCF and the --pass-only / \
             --min-af / --max-af / --min-depth filters if this is unexpected."
        );
    }

    // Fail early if no VCF contig matches any annotated gene — otherwise every
    // SNP lands outside every gene and pN/pS is NaN everywhere with no signal.
    {
        let gene_ids: HashSet<&str> = genes.iter().map(|g| g.seqid.as_str()).collect();
        let snp_chroms: HashSet<&str> = snps.iter().map(|s| s.chrom.as_str()).collect();
        if !snp_chroms.is_empty() && snp_chroms.is_disjoint(&gene_ids) {
            bail!(
                "No VCF CHROM matches any annotated gene sequence.\n  VCF uses e.g.:       {}\n  GFF uses:            {}\nContig names must match across reference, GFF, and VCF.",
                sample_names(&snp_chroms),
                sample_names(&gene_ids)
            );
        }
    }

    // Compute pN/pS
    let mut results = vcf_analysis::compute_pn_ps(
        &reference, &genes, &snps, gc, args.af_weighted, args.kappa, args.mk_fixed_af,
    );

    // Genome-wide pooled estimate and gene counts use ALL genes (pooling is
    // robust to low-count genes), captured before any --min-snps filtering.
    let total_genes = results.len();
    let genes_with_snps = results.iter().filter(|r| r.total_snps > 0.0).count();
    let genome_wide = vcf_analysis::genome_wide_pn_ps(&results);
    let gw_ci = if args.bootstrap > 0 {
        vcf_analysis::bootstrap_genome_wide_ci(&results, args.bootstrap, args.seed, 0.95)
    } else {
        None
    };

    // --min-snps drops unreliable low-count genes from the per-gene table,
    // plot, and neutrality test (but not from the pooled estimate above).
    let mut dropped = 0usize;
    if args.min_snps > 0 {
        let before = results.len();
        results.retain(|r| r.total_snps >= args.min_snps as f64);
        dropped = before - results.len();
        info!("--min-snps {}: kept {} of {} genes", args.min_snps, results.len(), before);
    }

    // Per-gene neutrality test correction (BH-FDR q-values + Bonferroni),
    // over the retained genes only.
    vcf_analysis::apply_multiple_testing(&mut results);

    // Write results
    let output_path = vcf_analysis::write_results(&results, &args.output, &args.format)?;
    info!("Results saved to {}", output_path);

    // McDonald-Kreitman test (optional).
    if args.mk {
        let mk_path = vcf_analysis::write_mk_results(&results, &args.output, &args.format)?;
        info!("McDonald-Kreitman results saved to {}", mk_path);
    }

    // Generate plots if requested
    if args.plot {
        let plot_path = vcf_analysis::write_pnps_plot(&results, &args.output, args.fdr)?;
        info!("Plot saved to {}", plot_path);
        let pv_path = vcf_analysis::write_pvalue_manhattan(&results, &args.output, args.fdr)?;
        info!("p-value Manhattan saved to {}", pv_path);
    }

    // Interactive HTML report.
    if args.report {
        let gc_label = format!("{} ({})", gc.id, gc.name);
        let rmeta = report::ReportMeta {
            n_samples,
            genetic_code: &gc_label,
            kappa: args.kappa,
            af_weighted: args.af_weighted,
            fdr: args.fdr,
            min_snps: args.min_snps,
            mk: args.mk,
            mk_fixed_af: args.mk_fixed_af,
            gw_ci,
        };
        let report_path =
            report::write_html_report(&results, genome_wide.as_ref(), &rmeta, &args.output)?;
        info!("Report saved to {}", report_path);
    }

    // Print summary statistics
    let (total_syn, total_nonsyn) = genome_wide
        .as_ref()
        .map(|gw| (gw.syn_snps, gw.nonsyn_snps))
        .unwrap_or((0.0, 0.0));
    let mode = if args.af_weighted { "πN/πS (AF-weighted)" } else { "pN/pS" };
    let ratio_name = if args.af_weighted { "πN/πS" } else { "pN/pS" };
    eprintln!("\n── {} Summary ──────────────────────────", mode);
    if args.kappa != 1.0 {
        eprintln!("  Site model:          ts/tv-weighted (kappa = {})", args.kappa);
    }
    eprintln!("  Genes analyzed:      {}", total_genes);
    eprintln!("  Genes with SNPs:     {}", genes_with_snps);
    if args.min_snps > 0 {
        eprintln!("  Genes kept (>= {} SNPs): {}  ({} dropped)", args.min_snps, results.len(), dropped);
    }
    eprintln!("  Total synonymous:    {:.2}", total_syn);
    eprintln!("  Total nonsynonymous: {:.2}", total_nonsyn);
    if let Some(gw) = &genome_wide {
        eprintln!("  ── Genome-wide (pooled) ──────────────────");
        eprintln!("  N / S sites:         {:.1} / {:.1}", gw.n_sites, gw.s_sites);
        eprintln!("  Overall pN / pS:     {:.6} / {:.6}", gw.pn, gw.ps);
        eprintln!("  Overall {}:       {}", ratio_name, vcf_analysis::format_ratio(gw.pn_ps));
        eprintln!("  Selection:           {}", vcf_analysis::selection_label(gw.pn_ps));
        if let Some((lo, hi)) = gw_ci {
            eprintln!("  95% CI ({} boot):    [{:.6}, {:.6}]", args.bootstrap, lo, hi);
        }
    }
    let n_tested = results.iter().filter(|r| r.p_value.is_finite()).count();
    if n_tested > 0 {
        let n_sig = results
            .iter()
            .filter(|r| r.q_value.is_finite() && r.q_value < args.fdr)
            .count();
        eprintln!("  ── Neutrality test (pN/pS = 1) ───────────");
        eprintln!("  Genes tested:        {}", n_tested);
        eprintln!("  Significant genes:   {}  (BH-FDR < {})", n_sig, args.fdr);
    } else if args.af_weighted {
        eprintln!("  (per-gene neutrality test skipped under --af-weighted)");
    }
    if args.mk {
        let with_data = results
            .iter()
            .filter(|r| r.mk_dn + r.mk_ds + r.mk_pn + r.mk_ps > 0)
            .count();
        let total_dn: u32 = results.iter().map(|r| r.mk_dn).sum();
        let total_ds: u32 = results.iter().map(|r| r.mk_ds).sum();
        let total_pn: u32 = results.iter().map(|r| r.mk_pn).sum();
        let total_ps: u32 = results.iter().map(|r| r.mk_ps).sum();
        eprintln!("  ── McDonald-Kreitman (fixed AF >= {}) ──", args.mk_fixed_af);
        eprintln!("  Genes with MK data:  {}", with_data);
        eprintln!("  Totals Dn/Ds/Pn/Ps:  {}/{}/{}/{}", total_dn, total_ds, total_pn, total_ps);
    }
    eprintln!("───────────────────────────────────────────");

    info!("VCF analysis completed successfully.");
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
