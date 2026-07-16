//! `eskaks vcf` subcommand: per-gene pN/pS orchestration.

use super::*;
use std::collections::{HashMap, HashSet};
use anyhow::{bail, Context};
use log::{info, warn};

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

/// Parse a per-gene divergence table (`gene<TAB or ,>dN/dS`, header optional)
/// for the report's polymorphism-vs-divergence panel.
fn parse_divergence(path: &str) -> anyhow::Result<HashMap<String, f64>> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("Failed to read --divergence file: {}", path))?;
    let mut map = HashMap::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split(['\t', ',']);
        let gene = match fields.next() {
            Some(g) => g.trim(),
            None => continue,
        };
        if let Some(v) = fields.next() {
            if let Ok(x) = v.trim().parse::<f64>() {
                map.insert(gene.to_string(), x);
            }
        }
    }
    if map.is_empty() {
        warn!("--divergence file {} yielded no gene→dN/dS pairs (expected 'gene<TAB>dN/dS')", path);
    } else {
        info!("Loaded {} per-gene divergence dN/dS values", map.len());
    }
    Ok(map)
}

pub(crate) fn run_vcf(args: cli::VcfArgs) -> anyhow::Result<()> {
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
    let genes_with_snps = results.iter().filter(|r| r.n_snps > 0).count();
    // Stratified core vs repetitive pooled ratios (always computed, for the
    // report). --exclude-repetitive makes the primary estimate core-only.
    let (core_gw, rep_gw) = vcf_analysis::genome_wide_core_repetitive(&results);
    let genome_wide = if args.exclude_repetitive {
        core_gw.clone()
    } else {
        vcf_analysis::genome_wide_pn_ps(&results)
    };
    let gw_ci = if args.bootstrap > 0 {
        // Resample the SAME gene set as the point estimate — core-only under
        // --exclude-repetitive — so the CI actually bounds the reported ratio.
        if args.exclude_repetitive {
            let core: Vec<_> = results.iter().filter(|r| !r.repetitive).cloned().collect();
            vcf_analysis::bootstrap_genome_wide_ci(&core, args.bootstrap, args.seed, 0.95)
        } else {
            vcf_analysis::bootstrap_genome_wide_ci(&results, args.bootstrap, args.seed, 0.95)
        }
    } else {
        None
    };

    // --min-snps drops unreliable low-count genes from the per-gene table,
    // plot, and neutrality test (but not from the pooled estimate above).
    let mut dropped = 0usize;
    if args.min_snps > 0 {
        let before = results.len();
        results.retain(|r| r.n_snps >= args.min_snps);
        dropped = before - results.len();
        info!("--min-snps {}: kept {} of {} genes", args.min_snps, results.len(), before);
    }

    // Per-gene neutrality test correction (BH-FDR q-values + Bonferroni),
    // over the retained genes only (excluding repetitive genes if requested).
    vcf_analysis::apply_multiple_testing(&mut results, args.exclude_repetitive);

    // Genomic-control diagnostic (always) and correction (opt-in). λ summarises
    // how far the tested p-values depart from the uniform null — inflated in
    // clonal, linked genomes.
    let lambda = vcf_analysis::genomic_inflation_lambda(&results);
    if args.genomic_control {
        vcf_analysis::apply_genomic_control(&mut results, lambda);
        info!("Genomic control applied (λ = {:.3})", lambda);
    }

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
        let command = std::env::args().collect::<Vec<_>>().join(" ");
        let vcf_file = args.vcf.join(", ");
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
            lambda,
            genomic_control: args.genomic_control,
            exclude_repetitive: args.exclude_repetitive,
            version: env!("CARGO_PKG_VERSION"),
            command: &command,
            vcf_file: &vcf_file,
            ref_file: &args.reference,
            gff_file: &args.gff,
            total_genes,
            genes_with_snps,
        };
        let divergence = match &args.divergence {
            Some(p) => Some(parse_divergence(p)?),
            None => None,
        };
        let report_path = report::write_html_report(
            &results, genome_wide.as_ref(), core_gw.as_ref(), rep_gw.as_ref(),
            &rmeta, divergence.as_ref(), &args.output,
        )?;
        info!("Report saved to {}", report_path);
    }
    if args.divergence.is_some() && !args.report {
        warn!("--divergence is only used by the interactive report; add --report or the file is ignored.");
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
    // Under --exclude-repetitive the pooled totals are core-only while the gene
    // counts above are over all genes; flag the scope so the block is self-consistent.
    let core_note = if args.exclude_repetitive { "   (core genes only)" } else { "" };
    eprintln!("  Total synonymous:    {:.2}{}", total_syn, core_note);
    eprintln!("  Total nonsynonymous: {:.2}{}", total_nonsyn, core_note);
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
    // "Tested" = the multiple-testing family (finite q), not merely genes with a raw
    // p — under --exclude-repetitive the repetitive genes are dropped from the family.
    let n_tested = results.iter().filter(|r| r.q_value.is_finite()).count();
    if n_tested > 0 {
        // Significance uses the GC-corrected q under --genomic-control, matching the report.
        let n_sig = results
            .iter()
            .filter(|r| {
                let q = if args.genomic_control { r.q_gc } else { r.q_value };
                q.is_finite() && q < args.fdr
            })
            .count();
        let corr = if args.genomic_control { "BH-FDR·GC" } else { "BH-FDR" };
        eprintln!("  ── Neutrality test (pN/pS = 1) ───────────");
        eprintln!("  Genes tested:        {}", n_tested);
        eprintln!("  Significant genes:   {}  ({} < {})", n_sig, corr, args.fdr);
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

