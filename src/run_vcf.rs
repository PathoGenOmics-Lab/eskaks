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
    // The AF filters are applied to the whole SNP set before the McDonald-Kreitman
    // fixed/polymorphic split, so they silently delete the variants MK depends on:
    // --max-af drops the fixed/divergence class (its help even recommends 0.99, which
    // equals the default --mk-fixed-af), and --min-af drops the polymorphic class.
    if args.mk && (args.min_af.is_some() || args.max_af.is_some()) {
        warn!(
            "--mk with --min-af/--max-af: the AF filter is applied BEFORE the \
             McDonald-Kreitman fixed/polymorphic split, so NI, alpha and Fisher's p are \
             computed from a truncated table (--max-af removes fixed/divergence variants, \
             --min-af removes polymorphisms). Interpret the MK columns with care."
        );
    }
    if (1..100).contains(&args.bootstrap) {
        warn!(
            "--bootstrap {} is very low; a 95% CI from so few replicates is unreliable (use >= 1000).",
            args.bootstrap
        );
    }

    crate::init_global_pool(args.workers);

    let ref_path = std::path::Path::new(&args.reference);
    let gff_path = std::path::Path::new(&args.gff);

    // The reference/GFF/VCF readers decompress gzip/bgzip transparently; guard
    // only against a directory path here for a clear early error.
    input::ensure_not_directory(&args.reference)?;
    input::ensure_not_directory(&args.gff)?;
    for v in &args.vcf {
        input::ensure_not_directory(v)?;
    }

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

    // n_effective is the number of sampled alleles for the diversity statistics:
    // the header sample-column count for a single multi-sample VCF, or the number
    // of files when merging one VCF per sample. 0 means "unknown" (AF-only input).
    let (snps, n_effective) = if n_samples == 1 {
        // Single VCF: use AF as-is (could be multi-sample VCF)
        let path = std::path::Path::new(&vcf_paths[0]);
        let snps = vcf::parse_vcf(path)?;
        info!("Found {} SNP records", snps.len());
        let n_eff = vcf::sample_count(path)?;
        (
            vcf::filter_snps(snps, args.pass_only, args.min_af, args.max_af, args.min_depth),
            n_eff,
        )
    } else {
        // Multiple single-sample VCFs: merge and compute AF as fraction of samples
        let merged = vcf::merge_vcfs(&vcf_paths, args.pass_only, args.min_depth)?;
        info!("Merged {} unique SNP positions from {} samples", merged.len(), n_samples);
        (
            vcf::filter_snps(merged, false, args.min_af, args.max_af, None),
            n_samples,
        )
    };
    info!("{} SNPs after filtering", snps.len());

    // What retaining per-sample carriers actually costs, reported rather than assumed:
    // one u64 per 64 samples per ALT allele. Absent (0 alleles) for an AF-only VCF with
    // no genotype columns, which is exactly the input that falls back to the AF bound.
    {
        let (mut alleles, mut bytes) = (0usize, 0usize);
        for cs in snps.iter().filter_map(|s| s.carriers.as_ref()) {
            alleles += cs.len();
            bytes += cs.iter().map(|c| c.heap_bytes()).sum::<usize>();
        }
        if alleles > 0 {
            info!(
                "Per-sample carriers retained for {} allele(s): {:.2} MiB ({} bytes)",
                alleles,
                bytes as f64 / (1024.0 * 1024.0),
                bytes
            );
        } else if n_effective > 0 {
            info!(
                "No per-sample carriers available; the same-codon check falls back to the \
                 allele-frequency bound."
            );
        }
    }

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
    let (mut results, diag): (Vec<_>, vcf_analysis::ComputeDiagnostics) = vcf_analysis::compute_pn_ps(
        &reference, &genes, &snps, gc, args.af_weighted, args.kappa, args.mk_fixed_af,
    );

    // ── Diagnostics: make an empty / garbage result impossible to mistake for
    // a clean run. ────────────────────────────────────────────────────────────
    // 1) VCF contigs that match no annotated gene: their SNPs were silently dropped.
    {
        let gene_ids: HashSet<&str> = genes.iter().map(|g| g.seqid.as_str()).collect();
        let mut unmatched: std::collections::BTreeMap<&str, usize> = std::collections::BTreeMap::new();
        for s in &snps {
            if !gene_ids.contains(s.chrom.as_str()) {
                *unmatched.entry(s.chrom.as_str()).or_insert(0) += 1;
            }
        }
        if !unmatched.is_empty() {
            let dropped: usize = unmatched.values().sum();
            let names: Vec<String> = unmatched.iter().map(|(c, n)| format!("{} ({} SNPs)", c, n)).collect();
            warn!(
                "{} of {} SNPs are on VCF contig(s) that match no annotated gene and were dropped: {}. \
                 Check that contig names agree across VCF, GFF and reference.",
                dropped, snps.len(), names.join(", ")
            );
        }
    }
    // 2) SNPs supplied but none landed inside a CDS: almost always a coordinate/build
    //    mismatch (contig names match but positions do not) rather than 'no data'.
    if !snps.is_empty() && diag.snps_in_cds == 0 {
        warn!(
            "{} SNPs were read but NONE fell inside an annotated CDS. Contig names match, so the \
             coordinates likely disagree with the GFF (wrong reference build, coordinate system, or \
             all variants are intergenic). Verify the VCF, GFF and reference share the same assembly.",
            snps.len()
        );
    }
    // 3) Reference genes that do not translate cleanly → wrong genetic code / frame.
    if diag.genes_with_internal_stops > 0 && !results.is_empty() {
        let frac = diag.genes_with_internal_stops as f64 / results.len().max(1) as f64;
        warn!(
            "{}/{} genes ({:.0}%) have an internal stop codon in their reference CDS — this usually \
             means the wrong --genetic-code (currently {}) or a frame/phase problem in the GFF.",
            diag.genes_with_internal_stops, results.len(), 100.0 * frac, args.genetic_code
        );
    }

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

    // Genome-wide (pooled) diversity, like the pooled pN/pS above, pools over ALL
    // genes, captured here BEFORE --min-snps drops low-count genes below, so the
    // headline π / θ_W / Tajima's D never silently depend on the --min-snps threshold.
    let genome_wide_div: Option<vcf_analysis::GenomeDiversity> = if args.diversity && n_effective >= 2
    {
        vcf_analysis::genome_wide_diversity(&results, n_effective)
    } else {
        None
    };

    // Per-codon recurrence scan (optional). Like the pooled estimates above, its
    // multiple-testing family is the whole coding genome, so it is built here from the
    // UNFILTERED results: --min-snps decides which genes are tabulated, not how many
    // codons a genome has. The parent-gene q-values are joined in after the gene-level
    // correction runs below.
    let mut codon_scan = if args.codon_scan {
        if n_effective > 0 && n_effective < 10 {
            warn!(
                "--codon-scan with only {} sample(s): recurrence is barely observable at this \
                 cohort size, since a residue has to collect several INDEPENDENT alleles before \
                 it can stand out. Expect an empty result rather than a negative one.",
                n_effective
            );
        }
        Some(vcf_analysis::compute_codon_scan(
            &results, gc, n_effective, args.exclude_repetitive,
        ))
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
        if results.is_empty() && before > 0 {
            warn!(
                "--min-snps {} dropped all {} genes — the per-gene table, plot and report will have \
                 no rows. Lower --min-snps. (The genome-wide pooled estimate still uses all genes.)",
                args.min_snps, before
            );
        }
    }

    // Per-gene neutrality test correction (BH-FDR q-values + Bonferroni),
    // over the retained genes only (excluding repetitive genes if requested).
    vcf_analysis::apply_multiple_testing(&mut results, args.exclude_repetitive);

    // The codon scan reports its parent gene's neutrality-test q-value as context, so
    // it can only be joined once that correction has been applied.
    if let Some(scan) = codon_scan.as_mut() {
        scan.attach_gene_qvalues(&results);
    }

    // Genomic-control diagnostic (always) and correction (opt-in). λ summarises
    // how far the tested p-values depart from the uniform null — inflated in
    // clonal, linked genomes.
    let lambda = vcf_analysis::genomic_inflation_lambda(&results);
    if args.genomic_control {
        vcf_analysis::apply_genomic_control(&mut results, lambda);
        info!("Genomic control applied (λ = {:.3})", lambda);
    }

    // Write results (collect every path so the summary can list them regardless of
    // log level — a user should always see where their output went).
    let mut written: Vec<String> = Vec::new();
    let output_path = vcf_analysis::write_results(&results, &args.output, &args.format)?;
    info!("Results saved to {}", output_path);
    written.push(output_path);

    // McDonald-Kreitman test (optional).
    if args.mk {
        let mk_path = vcf_analysis::write_mk_results(&results, &args.output, &args.format)?;
        info!("McDonald-Kreitman results saved to {}", mk_path);
        written.push(mk_path);
    }

    // Per-variant table (optional): the mutation-level detail behind each gene.
    if args.variants {
        let var_path =
            vcf_analysis::write_variants(&results, &args.output, &args.format, args.shared_codons)?;
        info!("Per-variant table saved to {}", var_path);
        written.push(var_path);
    }

    // Per-codon recurrence table (optional).
    if let Some(scan) = codon_scan.as_ref() {
        let codon_path = vcf_analysis::write_codon_scan(scan, &args.output, &args.format)?;
        info!("Per-codon recurrence table saved to {}", codon_path);
        written.push(codon_path);
    }

    // Population-diversity statistics (optional): π, Watterson θ, Tajima's D.
    if args.diversity {
        if n_effective >= 2 {
            let div_path =
                vcf_analysis::write_diversity(&results, n_effective, &args.output, &args.format)?;
            info!("Diversity table saved to {}", div_path);
            written.push(div_path);
            // Use the pre-filter genome-wide diversity captured above (pooled over ALL
            // genes), not a recompute over the --min-snps-filtered `results`.
            if let Some(gw) = genome_wide_div {
                // Explicitly requested via --diversity, so show the headline numbers.
                eprintln!("\n── Genome-wide diversity (n={}) ───────────", gw.n);
                eprintln!("  Segregating coding SNPs: {}", gw.s_seg);
                eprintln!("  piN / piS (per site):    {:.3e} / {:.3e}", gw.pi_n, gw.pi_s);
                eprintln!("  piN/piS:                 {:.4}", gw.pi_n_pi_s);
                eprintln!("  Watterson theta (site):  {:.3e}", gw.theta_w_per_site);
                eprintln!("  Tajima's D:              {:.4}", gw.tajima_d);
            }
        } else {
            warn!(
                "--diversity needs the sample size (a multi-sample VCF or several \
                 single-sample VCFs via --vcf-list); this input has no per-sample \
                 genotype columns, so π / θ_W / Tajima's D were skipped"
            );
        }
    }

    // Generate plots if requested. Each writer returns None when it has no data to plot
    // (e.g. no coding SNPs), so we only list files that were actually written.
    if args.plot {
        if let Some(plot_path) = vcf_analysis::write_pnps_plot(&results, &args.output, args.fdr)? {
            info!("Plot saved to {}", plot_path);
            written.push(plot_path);
        }
        if let Some(pv_path) = vcf_analysis::write_pvalue_manhattan(&results, &args.output, args.fdr)? {
            info!("p-value Manhattan saved to {}", pv_path);
            written.push(pv_path);
        }
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
            multi_snp_codons: diag.multi_snp_codons,
            cooccurring_codons: diag.cooccurring_codons,
            cooccurring_exact: diag.cooccurring_exact,
            mnv_alleles: diag.mnv_alleles,
            mnv_new_stops: diag.mnv_new_stops,
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
        written.push(report_path);
    }
    if args.divergence.is_some() && !args.report {
        warn!("--divergence is only used by the interactive report; add --report or the file is ignored.");
    }

    // Print the pN/pS summary: shown by default, hidden by --quiet, and forced by
    // --summary even under --quiet (symmetric with `eskaks fasta`).
    if args.summary || log::log_enabled!(log::Level::Warn) {
        let gw_totals = genome_wide.as_ref().map(|gw| (gw.syn_snps, gw.nonsyn_snps));
        let mode = if args.af_weighted { "πN/πS (AF-weighted)" } else { "pN/pS" };
        let ratio_name = if args.af_weighted { "πN/πS" } else { "pN/pS" };
        eprintln!("\n── {} Summary ──────────────────────────", mode);
        if args.kappa != 1.0 {
            eprintln!("  Site model:          ts/tv-weighted (kappa = {})", args.kappa);
        }
        eprintln!("  Genes analyzed:      {}", total_genes);
        eprintln!("  Genes with SNPs:     {}", genes_with_snps);
        eprintln!("  SNPs used (in CDS):  {} of {} parsed", diag.snps_in_cds, snps.len());
        if diag.ref_mismatch > 0 {
            eprintln!("  SNPs skipped (REF≠ref): {}", diag.ref_mismatch);
        }
        // Codons carrying more than one SNP. The indented subset is the one where the SNPs
        // provably occur together in one genome. Where per-sample genotypes survived
        // parsing that subset is exact AND the codon those samples carry is what was
        // scored; otherwise it rests on the allele-frequency bound, which is a floor and
        // leaves the joint change unevaluated. The basis is named rather than left to
        // guess, because the two lead to different numbers.
        if diag.multi_snp_codons > 0 {
            eprintln!("  Codons with >1 SNP:  {}", diag.multi_snp_codons);
            let basis = if diag.codons_with_carriers == diag.multi_snp_codons {
                "exact, from per-sample genotypes"
            } else if diag.codons_with_carriers == 0 {
                "allele-frequency bound: a floor, not a count"
            } else {
                "part exact, part allele-frequency bound"
            };
            eprintln!("    checked:           {}", basis);
            if diag.cooccurring_codons > 0 {
                eprintln!(
                    "    shared by a sample: {}  (in {} gene(s))",
                    diag.cooccurring_codons, diag.genes_with_cooccurring
                );
                let bound = diag.cooccurring_codons - diag.cooccurring_exact;
                if diag.cooccurring_exact > 0 && bound > 0 {
                    eprintln!(
                        "      observed / bound: {} / {}",
                        diag.cooccurring_exact, bound
                    );
                }
                if diag.mnv_alleles > 0 {
                    eprintln!(
                        "      scored jointly:   {} allele(s), by Nei-Gojobori pathway averaging",
                        diag.mnv_alleles
                    );
                }
                if diag.mnv_excluded_from_mk > 0 {
                    eprintln!(
                        "      out of MK / SFS:  {}  (both need whole alleles)",
                        diag.mnv_excluded_from_mk
                    );
                }
                if diag.mnv_new_stops > 0 {
                    eprintln!(
                        "      premature stops:  {} codon(s) stop only when taken together",
                        diag.mnv_new_stops
                    );
                }
                if bound > 0 {
                    eprintln!(
                        "      not evaluated:    {} codon(s): no genotypes, so no codon to score",
                        bound
                    );
                }
            }
        }
        if args.min_snps > 0 {
            eprintln!("  Genes kept (>= {} SNPs): {}  ({} dropped)", args.min_snps, results.len(), dropped);
        }
        // Under --exclude-repetitive the pooled totals are core-only while the gene
        // counts above are over all genes; flag the scope so the block is self-consistent.
        let core_note = if args.exclude_repetitive { "   (core genes only)" } else { "" };
        match gw_totals {
            Some((total_syn, total_nonsyn)) => {
                eprintln!("  Total synonymous:    {:.2}{}", total_syn, core_note);
                eprintln!("  Total nonsynonymous: {:.2}{}", total_nonsyn, core_note);
            }
            None => {
                // No pooled estimate: report n/a for the totals, consistent with the pooled
                // ratio line below, rather than a misleading 0.00 that reads as "zero coding
                // SNPs" when in fact every contributing gene was filtered out.
                let why = if args.exclude_repetitive {
                    "   (all genes excluded as repetitive)"
                } else {
                    "   (no gene contributed a coding SNP)"
                };
                eprintln!("  Total synonymous:    n/a{}", why);
                eprintln!("  Total nonsynonymous: n/a{}", why);
            }
        }
        if let Some(gw) = &genome_wide {
            eprintln!("  ── Genome-wide (pooled) ──────────────────");
            eprintln!("  N / S sites:         {:.1} / {:.1}", gw.n_sites, gw.s_sites);
            eprintln!("  Overall pN / pS:     {:.6} / {:.6}", gw.pn, gw.ps);
            eprintln!("  Overall {}:       {}", ratio_name, vcf_analysis::format_ratio(gw.pn_ps));
            eprintln!("  Selection:           {}", vcf_analysis::selection_label(gw.pn_ps));
            if let Some((lo, hi)) = gw_ci {
                eprintln!("  95% CI ({} boot):    [{:.6}, {:.6}]", args.bootstrap, lo, hi);
            }
        } else {
            // No pooled estimate: say so explicitly rather than leaving a silent gap.
            eprintln!("  ── Genome-wide (pooled) ──────────────────");
            let why = if args.exclude_repetitive {
                " (all genes were excluded as repetitive — drop --exclude-repetitive?)"
            } else {
                " (no gene contributed any coding SNP)"
            };
            eprintln!("  Overall {}:       n/a{}", ratio_name, why);
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
        if let Some(scan) = codon_scan.as_ref() {
            // m and theta are the whole null: without them the p-values cannot be
            // reproduced or sanity-checked, so they belong in the summary, not only
            // in the docs.
            eprintln!("  ── Per-codon recurrence scan ─────────────");
            eprintln!("  Codons with SNPs:    {}", scan.rows.len());
            eprintln!(
                "  Family (m):          {} codons, {} possible nonsyn changes{}",
                scan.family_m, scan.family_poss_nonsyn,
                if args.exclude_repetitive { "  (core genes only)" } else { "" }
            );
            eprintln!(
                "  Distinct nonsyn alleles: {}  (theta = {:.3e} per possible change)",
                scan.observed_nonsyn_alleles, scan.theta
            );
            if scan.n_samples == 0 {
                // The test only needs allele identity, so it still ran; only the
                // descriptive carrier columns are unavailable.
                eprintln!("  Carrier columns:     NA (input has no genotypes or sample count)");
            }
            if scan.suppressed_cooccurring > 0 {
                eprintln!(
                    "  Not tested (one haplotype): {}",
                    scan.suppressed_cooccurring
                );
            }
            let sig = scan.significant(args.fdr);
            eprintln!("  Significant codons:  {}  (BH-FDR < {})", sig, args.fdr);
            if sig == 0 {
                eprintln!("    (no residue collected more independent alleles than chance allows)");
            }
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
        // Tell the user where the output went (info-level "saved to" lines are hidden at
        // the default log level). One file per line, matching the fasta "Done" block.
        if written.is_empty() {
            eprintln!("  Output:     (none)");
        } else {
            eprintln!("  Output:");
            for p in &written {
                eprintln!("    {}", p);
            }
        }
        eprintln!("───────────────────────────────────────────");
    }

    info!("VCF analysis completed successfully.");
    Ok(())
}

