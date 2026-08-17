//! Per-gene pN/pS and McDonald-Kreitman table writers.

use super::*;
// Gene/chrom names come from GFF3 attributes (percent-decoded), so they can carry
// quotes, delimiters, or control chars; escape/quote them like every other writer.
use crate::textfmt::{name_field, json_escape};

/// Write pN/pS results to a file.
pub fn write_results(
    results: &[GenePnPs],
    prefix: &str,
    format: &crate::models::OutputFormat,
) -> anyhow::Result<String> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let ext = format.extension();
    let output_path = format!("{}_pnps.{}", prefix, ext);

    match format {
        crate::models::OutputFormat::Json => {
            let mut file = BufWriter::new(File::create(&output_path)?);
            writeln!(file, "[")?;
            for (i, r) in results.iter().enumerate() {
                let comma = if i + 1 < results.len() { "," } else { "" };
                writeln!(
                    file,
                    "  {{\"gene\":\"{}\",\"chrom\":\"{}\",\"start\":{},\"end\":{},\"strand\":\"{}\",\"length_bp\":{},\"N_sites\":{:.4},\"S_sites\":{:.4},\"exp_N_frac\":{},\"pN\":{},\"pS\":{},\"pN_pS\":{},\"pN_pS_lo\":{},\"pN_pS_hi\":{},\"nonsyn_snps\":{:.4},\"syn_snps\":{:.4},\"total_snps\":{:.4},\"p_value\":{},\"q_value_bh\":{},\"p_bonferroni\":{},\"p_gc\":{},\"q_gc_bh\":{}}}{}",
                    json_escape(&r.name), json_escape(&r.chrom), r.genome_start, r.genome_end, r.strand, r.length_bp,
                    r.n_sites, r.s_sites, format_json_num(exp_n_frac(r)),
                    format_json_f64(r.pn), format_json_f64(r.ps), format_json_f64(r.pn_ps),
                    format_json_num(r.pn_ps_lo), format_json_num(r.pn_ps_hi),
                    r.nonsyn_snps, r.syn_snps, r.total_snps,
                    format_json_num(r.p_value), format_json_num(r.q_value), format_json_num(r.p_bonferroni),
                    format_json_num(r.p_gc), format_json_num(r.q_gc),
                    comma
                )?;
            }
            writeln!(file, "]")?;
        }
        _ => {
            let sep = format.separator();
            let mut file = BufWriter::new(File::create(&output_path)?);
            // New columns (Wilson CI + genomic-control p/q) are appended at the end so
            // existing 0-indexed column positions stay stable for downstream parsers.
            writeln!(
                file,
                "Gene{s}Length_bp{s}N_sites{s}S_sites{s}pN{s}pS{s}pN/pS{s}Nonsyn_SNPs{s}Syn_SNPs{s}Total_SNPs{s}Chrom{s}Start{s}End{s}Strand{s}Exp_N_frac{s}P_value{s}Q_value_BH{s}P_Bonferroni{s}pN/pS_lo{s}pN/pS_hi{s}P_GC{s}Q_GC_BH",
                s = sep
            )?;
            for r in results {
                writeln!(
                    file,
                    "{}{s}{}{s}{:.4}{s}{:.4}{s}{:.6}{s}{:.6}{s}{}{s}{:.4}{s}{:.4}{s}{:.4}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}",
                    name_field(&r.name, sep), r.length_bp, r.n_sites, r.s_sites,
                    r.pn, r.ps, format_ratio(r.pn_ps),
                    r.nonsyn_snps, r.syn_snps, r.total_snps,
                    name_field(&r.chrom, sep), r.genome_start, r.genome_end, r.strand,
                    format_pval(exp_n_frac(r)), format_pval(r.p_value),
                    format_pval(r.q_value), format_pval(r.p_bonferroni),
                    format_ratio(r.pn_ps_lo), format_ratio(r.pn_ps_hi),
                    format_pval(r.p_gc), format_pval(r.q_gc),
                    s = sep
                )?;
            }
        }
    }

    Ok(output_path)
}

/// Write the per-coding-SNP variants table behind the pN/pS counts: one row per
/// ALT allele located in a CDS, with its genomic position, base change, protein
/// residue and amino-acid change (e.g. `S315T`), allele frequency, and effect
/// (synonymous / missense / nonsense / stop_loss). This is the mutation-level key
/// a user joins to the WHO catalogue or TB-Profiler; nonsense/stop-loss changes
/// are included even though they are excluded from the pN/pS site & SNP counts.
pub fn write_variants(
    results: &[GenePnPs],
    prefix: &str,
    format: &crate::models::OutputFormat,
) -> anyhow::Result<String> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let ext = format.extension();
    let output_path = format!("{}_variants.{}", prefix, ext);
    let mut file = BufWriter::new(File::create(&output_path)?);
    let change = |v: &Variant| format!("{}{}{}", v.ref_aa as char, v.aa_pos, v.alt_aa as char);

    match format {
        crate::models::OutputFormat::Json => {
            writeln!(file, "[")?;
            let mut first = true;
            for g in results {
                for v in &g.variants {
                    let comma = if first { "" } else { ",\n" };
                    first = false;
                    write!(
                        file,
                        "{comma}  {{\"gene\":\"{}\",\"chrom\":\"{}\",\"pos\":{},\"strand\":\"{}\",\"ref\":\"{}\",\"alt\":\"{}\",\"aa_pos\":{},\"ref_aa\":\"{}\",\"alt_aa\":\"{}\",\"change\":\"{}\",\"af\":{},\"effect\":\"{}\"}}",
                        json_escape(&g.name), json_escape(&g.chrom), v.pos, g.strand,
                        v.ref_allele as char, v.alt_allele as char, v.aa_pos,
                        v.ref_aa as char, v.alt_aa as char, change(v),
                        format_json_num(v.af), v.effect.label()
                    )?;
                }
            }
            writeln!(file, "\n]")?;
        }
        _ => {
            let sep = format.separator();
            writeln!(
                file,
                "Gene{s}Chrom{s}Pos{s}Strand{s}Ref{s}Alt{s}AA_Pos{s}Ref_AA{s}Alt_AA{s}Change{s}AF{s}Effect",
                s = sep
            )?;
            for g in results {
                for v in &g.variants {
                    writeln!(
                        file,
                        "{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{:.4}{s}{}",
                        name_field(&g.name, sep), name_field(&g.chrom, sep), v.pos, g.strand,
                        v.ref_allele as char, v.alt_allele as char, v.aa_pos,
                        v.ref_aa as char, v.alt_aa as char, change(v), v.af, v.effect.label(),
                        s = sep
                    )?;
                }
            }
        }
    }

    Ok(output_path)
}

/// Write the per-codon recurrence table (`<output>_codons.<ext>`): one row per codon
/// carrying at least one coding SNP, ranked by how unlikely its number of **distinct**
/// nonsynonymous alleles is under a genome-wide per-change rate.
///
/// The only tested column is `P_Recurrence` (with its BH q-value). Carrier counts and
/// allele frequencies are descriptive and deliberately carry no p-value: without a
/// phylogeny, one old mutation in an expanded lineage cannot be told from many
/// independent origins. `Gene_pN_pS` and `Gene_Q_BH` are joined from the parent gene so
/// a globally elevated gene is not misread as a codon-specific hit. See
/// [`crate::vcf_analysis::compute_codon_scan`] for the null and its limits.
pub fn write_codon_scan(
    scan: &CodonScan,
    prefix: &str,
    format: &crate::models::OutputFormat,
) -> anyhow::Result<String> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let ext = format.extension();
    let output_path = format!("{}_codons.{}", prefix, ext);
    let mut file = BufWriter::new(File::create(&output_path)?);

    // "NA" for a codon with no amino-acid change, so an empty cell is never ambiguous.
    fn changes(r: &CodonRow) -> &str {
        if r.aa_changes.is_empty() {
            "NA"
        } else {
            r.aa_changes.as_str()
        }
    }
    let codon = |r: &CodonRow| String::from_utf8_lossy(&r.ref_codon).into_owned();

    match format {
        crate::models::OutputFormat::Json => {
            writeln!(file, "[")?;
            for (i, r) in scan.rows.iter().enumerate() {
                let comma = if i + 1 < scan.rows.len() { "," } else { "" };
                let carriers = |v: Option<usize>| match v {
                    Some(c) => c.to_string(),
                    None => "null".to_string(),
                };
                // A codon with no amino-acid change is `null`, not the string "NA":
                // JSON's rule everywhere else in eskaks is that an absent value is null.
                let aa_changes = if r.aa_changes.is_empty() {
                    "null".to_string()
                } else {
                    format!("\"{}\"", json_escape(&r.aa_changes))
                };
                writeln!(
                    file,
                    "  {{\"gene\":\"{}\",\"chrom\":\"{}\",\"strand\":\"{}\",\"aa_pos\":{},\"ref_aa\":\"{}\",\"ref_codon\":\"{}\",\"poss_nonsyn\":{},\"poss_syn\":{},\"nonsyn_alleles\":{},\"syn_alleles\":{},\"nonsense_alleles\":{},\"distinct_aa\":{},\"aa_changes\":{},\"carriers_max\":{},\"carriers_sum\":{},\"max_af\":{},\"exp_nonsyn_alleles\":{},\"p_recurrence\":{},\"q_recurrence_bh\":{},\"cooccurring\":{},\"repetitive\":{},\"gene_pn_ps\":{},\"gene_q_bh\":{}}}{}",
                    json_escape(&r.gene), json_escape(&r.chrom), r.strand, r.aa_pos,
                    r.ref_aa as char, codon(r),
                    r.poss_nonsyn, r.poss_syn,
                    r.nonsyn_alleles, r.syn_alleles, r.nonsense_alleles,
                    r.distinct_aa, aa_changes,
                    carriers(r.carriers_max), carriers(r.carriers_sum),
                    format_json_num(r.max_af), format_json_num(r.exp_nonsyn_alleles),
                    format_json_num(r.p_recurrence), format_json_num(r.q_recurrence),
                    r.cooccurring, r.repetitive,
                    format_json_num(r.gene_pn_ps), format_json_num(r.gene_q_bh),
                    comma
                )?;
            }
            writeln!(file, "]")?;
        }
        _ => {
            let sep = format.separator();
            writeln!(
                file,
                "Gene{s}Chrom{s}Strand{s}AA_Pos{s}Ref_AA{s}Ref_Codon{s}Poss_Nonsyn{s}Poss_Syn{s}Nonsyn_Alleles{s}Syn_Alleles{s}Nonsense_Alleles{s}Distinct_AA{s}AA_Changes{s}Carriers_Max{s}Carriers_Sum{s}Max_AF{s}Exp_Nonsyn_Alleles{s}P_Recurrence{s}Q_Recurrence_BH{s}Cooccurring{s}Repetitive{s}Gene_pN_pS{s}Gene_Q_BH",
                s = sep
            )?;
            let carriers = |v: Option<usize>| match v {
                Some(c) => c.to_string(),
                None => "NA".to_string(),
            };
            // 4 dp like the variants table's AF, but "NA" rather than a bare "NaN" when
            // no allele carried a usable frequency.
            let af_cell = |v: f64| {
                if v.is_nan() {
                    "NA".to_string()
                } else {
                    format!("{:.4}", crate::textfmt::norm_zero(v))
                }
            };
            for r in &scan.rows {
                writeln!(
                    file,
                    "{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}",
                    name_field(&r.gene, sep), name_field(&r.chrom, sep), r.strand, r.aa_pos,
                    r.ref_aa as char, codon(r),
                    r.poss_nonsyn, r.poss_syn,
                    r.nonsyn_alleles, r.syn_alleles, r.nonsense_alleles,
                    r.distinct_aa, name_field(changes(r), sep),
                    carriers(r.carriers_max), carriers(r.carriers_sum),
                    af_cell(r.max_af),
                    format_pval(r.exp_nonsyn_alleles),
                    format_pval(r.p_recurrence), format_pval(r.q_recurrence),
                    r.cooccurring, r.repetitive,
                    format_ratio(r.gene_pn_ps), format_pval(r.gene_q_bh),
                    s = sep
                )?;
            }
        }
    }

    Ok(output_path)
}

/// Keep a derived count only if the site is genuinely *segregating*: `None` when it is
/// monomorphic (`k == 0` absent, or `k >= n` fixed for the ALT), so it is excluded
/// from π / θ_W / Tajima's D, never clamped into range.
fn segregating(k: usize, n: usize) -> Option<usize> {
    if k == 0 || k >= n {
        None
    } else {
        Some(k)
    }
}

/// A segregating site: its derived-allele count paired with the sample size actually
/// observed there (the per-site called count, or the nominal n for an AF-only record).
type SegSite = (usize, usize);

/// Segregating derived-allele counts of a variant set, each paired with the sample
/// size actually observed at its site (per-site called count from the GT columns, or
/// the nominal `n` for an AF-only / merged record). Split into (synonymous, missense);
/// monomorphic sites (absent or fixed for the ALT) and nonsense/stop-loss are dropped.
/// The per-site denominator is what makes a no-call-heavy site score as fixed or
/// segregating consistently with the AF path, instead of diluting its frequency with
/// samples that were never genotyped.
fn segregating_counts(
    variants: &[Variant],
    n: usize,
) -> (Vec<SegSite>, Vec<SegSite>) {
    let (mut syn, mut mis) = (Vec::new(), Vec::new());
    for v in variants {
        let bucket = match v.effect {
            SnpEffect::Synonymous => &mut syn,
            SnpEffect::Missense => &mut mis,
            _ => continue, // nonsense/stop-loss excluded, as in the pN/pS counts
        };
        // Exact GT-derived count over its own called sample size; fall back to
        // round(af * n) over the nominal n only when the record had no genotypes (a
        // merged multi-VCF, where af = count / n is exact, or an AF-only VCF).
        let (derived, ni) = match v.gt_derived {
            Some(g) => (g, v.gt_called.unwrap_or(n)),
            None => ((v.af * n as f64).round() as usize, n),
        };
        if let Some(k) = segregating(derived, ni) {
            bucket.push((k, ni));
        }
    }
    (syn, mis)
}

/// The most common per-site called sample size among the given sites, ties broken
/// toward the larger n. `None` for an empty set. Used as the single representative
/// sample size for θ_W and Tajima's D, which (unlike π) are not defined per-site.
fn modal_called<'a>(sites: impl Iterator<Item = &'a SegSite>) -> Option<usize> {
    let mut freq: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for &(_, ni) in sites {
        *freq.entry(ni).or_insert(0) += 1;
    }
    freq.into_iter()
        .max_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)))
        .map(|(ni, _)| ni)
}

/// A computed diversity summary from split segregating counts and site totals.
struct DiversityRow {
    s_seg: usize,
    pi_n: f64,
    pi_s: f64,
    pi_ratio: f64,
    theta_w: f64,
    tajima_d: f64,
}

fn diversity_row(
    syn: &[SegSite],
    mis: &[SegSite],
    n_sites: f64,
    s_sites: f64,
    n: usize,
) -> DiversityRow {
    let s_seg = syn.len() + mis.len();
    // π keeps each site's own sample size n_i. θ_W and Tajima's D take a single sample
    // size, so use the most common per-site called count (n_eff), falling back to the
    // nominal n when nothing segregates. With no missing data every n_i == n, so this
    // reduces exactly to the fixed-n formulas.
    let n_eff = modal_called(syn.iter().chain(mis.iter())).unwrap_or(n);
    let pi_n = if n_sites > 0.0 { crate::stats::theta_pi_varn(mis) / n_sites } else { f64::NAN };
    let pi_s = if s_sites > 0.0 { crate::stats::theta_pi_varn(syn) / s_sites } else { f64::NAN };
    // Mirror the pN/pS ratio guard (see pnps.rs): a positive πN over zero πS is a
    // divergent ratio (+∞, all nonsynonymous diversity), not "no data". Collapsing it to
    // NaN would hide that signal and disagree with the pN/pS table on the same situation.
    // Only πN == πS == 0 (or an absent site denominator) is genuinely undefined.
    let pi_ratio = if pi_s > 0.0 {
        pi_n / pi_s
    } else if pi_n > 0.0 {
        f64::INFINITY
    } else {
        f64::NAN
    };
    let total = n_sites + s_sites;
    let theta_w =
        if total > 0.0 { crate::stats::theta_watterson(n_eff, s_seg) / total } else { f64::NAN };
    let all: Vec<SegSite> = syn.iter().chain(mis.iter()).copied().collect();
    let tajima_d = crate::stats::tajimas_d(n_eff, s_seg, crate::stats::theta_pi_varn(&all));
    DiversityRow { s_seg, pi_n, pi_s, pi_ratio, theta_w, tajima_d }
}

/// Genome-wide (pooled) diversity summary over all coding SNPs.
#[derive(Debug, Clone, Copy)]
pub struct GenomeDiversity {
    pub n: usize,
    pub s_seg: usize,
    pub pi_n: f64,
    pub pi_s: f64,
    pub pi_n_pi_s: f64,
    pub theta_w_per_site: f64,
    pub tajima_d: f64,
}

/// Pool every gene's coding SNPs into a single genome-wide diversity summary:
/// per-site πN and πS (and their ratio), per-site Watterson θ, and Tajima's D
/// over all coding sites. Needs the sample size `n ≥ 2` (returns None otherwise).
pub fn genome_wide_diversity(results: &[GenePnPs], n: usize) -> Option<GenomeDiversity> {
    if n < 2 {
        return None;
    }
    let (mut syn, mut mis) = (Vec::new(), Vec::new());
    let (mut n_sites, mut s_sites) = (0.0f64, 0.0f64);
    for g in results {
        n_sites += g.n_sites;
        s_sites += g.s_sites;
        let (gs, gm) = segregating_counts(&g.variants, n);
        syn.extend(gs);
        mis.extend(gm);
    }
    let r = diversity_row(&syn, &mis, n_sites, s_sites, n);
    Some(GenomeDiversity {
        n,
        s_seg: r.s_seg,
        pi_n: r.pi_n,
        pi_s: r.pi_s,
        pi_n_pi_s: r.pi_ratio,
        theta_w_per_site: r.theta_w,
        tajima_d: r.tajima_d,
    })
}

/// Write the per-gene diversity table (`<output>_diversity.<ext>`): sample size,
/// segregating coding SNPs, per-site πN and πS (nucleotide diversity, the
/// within-species analogue of pN/pS) and their ratio, per-site Watterson θ, and
/// Tajima's D — the SFS neutrality test (D < 0: excess of rare variants; D > 0:
/// intermediate-frequency excess). Requires the sample size `n ≥ 2`, and assumes
/// haploid genotypes (as for *M. tuberculosis*): a site's derived count is taken
/// straight from the GT columns when present (exact, and independent of any INFO/AF
/// that disagrees with them), falling back to round(AF·n) only for a merged multi-VCF
/// (where AF = carriers / n is itself exact) or an AF-only VCF. Multiallelic sites
/// contribute one entry per ALT rather than a single site-spectrum term, so π and the
/// segregating-site count are mildly inflated at such sites (rare in clonal genomes);
/// overlapping genes count a shared SNP once per gene.
pub fn write_diversity(
    results: &[GenePnPs],
    n: usize,
    prefix: &str,
    format: &crate::models::OutputFormat,
) -> anyhow::Result<String> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let ext = format.extension();
    let output_path = format!("{}_diversity.{}", prefix, ext);
    let mut file = BufWriter::new(File::create(&output_path)?);

    match format {
        crate::models::OutputFormat::Json => {
            writeln!(file, "[")?;
            for (i, g) in results.iter().enumerate() {
                let (syn, mis) = segregating_counts(&g.variants, n);
                let d = diversity_row(&syn, &mis, g.n_sites, g.s_sites, n);
                let comma = if i + 1 < results.len() { "," } else { "" };
                writeln!(
                    file,
                    "  {{\"gene\":\"{}\",\"chrom\":\"{}\",\"n_samples\":{},\"s_seg\":{},\"piN\":{},\"piS\":{},\"piN_piS\":{},\"theta_w\":{},\"tajima_d\":{}}}{}",
                    json_escape(&g.name), json_escape(&g.chrom), n, d.s_seg,
                    format_json_num(d.pi_n), format_json_num(d.pi_s), format_json_num(d.pi_ratio),
                    format_json_num(d.theta_w), format_json_num(d.tajima_d), comma
                )?;
            }
            writeln!(file, "]")?;
        }
        _ => {
            let sep = format.separator();
            writeln!(
                file,
                "Gene{s}Chrom{s}N_samples{s}S_seg{s}piN{s}piS{s}piN/piS{s}Theta_W{s}Tajima_D",
                s = sep
            )?;
            for g in results {
                let (syn, mis) = segregating_counts(&g.variants, n);
                let d = diversity_row(&syn, &mis, g.n_sites, g.s_sites, n);
                writeln!(
                    file,
                    "{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}",
                    name_field(&g.name, sep), name_field(&g.chrom, sep), n, d.s_seg,
                    format_pval(d.pi_n), format_pval(d.pi_s), format_ratio(d.pi_ratio),
                    format_pval(d.theta_w), format_ratio(d.tajima_d),
                    s = sep
                )?;
            }
        }
    }
    Ok(output_path)
}

/// Write per-gene McDonald-Kreitman results: the 2×2 fixed/polymorphic table
/// (Dn, Ds, Pn, Ps), Neutrality Index, alpha (proportion of adaptive
/// substitutions), a two-sided Fisher exact p-value, and its BH-FDR q-value.
///
/// "Fixed" means AF >= the chosen threshold *within the sample* (reference-
/// polarized MK); this conflates high-frequency derived alleles with true
/// between-species divergence, so it is a screen, not a substitute for an
/// outgroup-based MK test.
pub fn write_mk_results(
    results: &[GenePnPs],
    prefix: &str,
    format: &crate::models::OutputFormat,
) -> anyhow::Result<String> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    // Only genes with at least one classified difference are informative.
    let genes: Vec<&GenePnPs> = results
        .iter()
        .filter(|r| r.mk_dn + r.mk_ds + r.mk_pn + r.mk_ps > 0)
        .collect();

    // Two-sided Fisher p per gene, then Benjamini-Hochberg across tested genes.
    let pvals: Vec<f64> = genes
        .iter()
        .map(|r| {
            crate::stats::fisher_exact_two_sided(
                r.mk_dn as u64,
                r.mk_ds as u64,
                r.mk_pn as u64,
                r.mk_ps as u64,
            )
        })
        .collect();
    let qvals = crate::stats::benjamini_hochberg(&pvals);

    let ext = format.extension();
    let output_path = format!("{}_mk.{}", prefix, ext);
    let mut file = BufWriter::new(File::create(&output_path)?);

    let ni = |r: &GenePnPs| {
        let (dn, ds, pn, ps) = (r.mk_dn as f64, r.mk_ds as f64, r.mk_pn as f64, r.mk_ps as f64);
        if ps > 0.0 && dn > 0.0 {
            (pn * ds) / (ps * dn)
        } else {
            f64::NAN
        }
    };
    let alpha = |r: &GenePnPs| {
        let (dn, ds, pn, ps) = (r.mk_dn as f64, r.mk_ds as f64, r.mk_pn as f64, r.mk_ps as f64);
        if dn > 0.0 && ps > 0.0 {
            1.0 - (ds * pn) / (dn * ps)
        } else {
            f64::NAN
        }
    };

    if let crate::models::OutputFormat::Json = format {
        writeln!(file, "[")?;
        for (i, r) in genes.iter().enumerate() {
            let comma = if i + 1 < genes.len() { "," } else { "" };
            writeln!(
                file,
                "  {{\"gene\":\"{}\",\"chrom\":\"{}\",\"start\":{},\"end\":{},\"strand\":\"{}\",\"Dn\":{},\"Ds\":{},\"Pn\":{},\"Ps\":{},\"NI\":{},\"alpha\":{},\"fisher_p\":{},\"fisher_q_bh\":{}}}{}",
                json_escape(&r.name), json_escape(&r.chrom), r.genome_start, r.genome_end, r.strand,
                r.mk_dn, r.mk_ds, r.mk_pn, r.mk_ps,
                format_json_num(ni(r)), format_json_num(alpha(r)),
                format_json_num(pvals[i]), format_json_num(qvals[i]), comma
            )?;
        }
        writeln!(file, "]")?;
    } else {
        let sep = format.separator();
        writeln!(
            file,
            "Gene{s}Chrom{s}Start{s}End{s}Strand{s}Dn{s}Ds{s}Pn{s}Ps{s}NI{s}alpha{s}Fisher_p{s}Fisher_q_BH",
            s = sep
        )?;
        for (i, r) in genes.iter().enumerate() {
            writeln!(
                file,
                "{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}",
                name_field(&r.name, sep), name_field(&r.chrom, sep), r.genome_start, r.genome_end, r.strand,
                r.mk_dn, r.mk_ds, r.mk_pn, r.mk_ps,
                format_pval(ni(r)), format_pval(alpha(r)),
                format_pval(pvals[i]), format_pval(qvals[i]),
                s = sep
            )?;
        }
    }

    Ok(output_path)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn var(effect: SnpEffect, gt_derived: Option<usize>, gt_called: Option<usize>, af: f64) -> Variant {
        Variant {
            pos: 1, ref_allele: b'A', alt_allele: b'C', aa_pos: 1, ref_codon: *b"GCT",
            ref_aa: b'A', alt_aa: b'C', af, gt_derived, gt_called, effect,
        }
    }

    #[test]
    fn diversity_scores_no_calls_on_their_called_sample() {
        // 2 derived of 2 CALLED (the other 2 of a 4-sample cohort are no-calls) is fixed
        // among the genotyped samples, exactly as its AF = 1.0 says, so it must NOT count
        // as segregating even though the nominal cohort is n = 4. This is Bug B.
        let (syn, mis) = segregating_counts(&[var(SnpEffect::Synonymous, Some(2), Some(2), 1.0)], 4);
        assert!(syn.is_empty() && mis.is_empty(), "no-call-fixed site must be dropped");

        // Its genuinely-called counterpart (2 derived of 4 called) IS segregating, and
        // carries its own per-site sample size.
        let (syn, _) = segregating_counts(&[var(SnpEffect::Synonymous, Some(2), Some(4), 0.5)], 4);
        assert_eq!(syn, vec![(2usize, 4usize)]);

        // No genotype columns (AF-only / merged): fall back to round(af * nominal n).
        let (_, mis) = segregating_counts(&[var(SnpEffect::Missense, None, None, 0.5)], 4);
        assert_eq!(mis, vec![(2usize, 4usize)]);
    }

    #[test]
    fn pi_ratio_diverges_to_infinity_when_only_nonsynonymous_diversity() {
        // πN > 0 with πS == 0 (synonymous sites exist but none segregate) is a divergent
        // ratio, not "no data"; it must mirror the pN/pS table, which reports +inf here.
        let seg = vec![(2usize, 4usize)];
        let r = diversity_row(&[], &seg, 10.0, 10.0, 4);
        assert!(r.pi_n > 0.0 && r.pi_s == 0.0);
        assert!(r.pi_ratio.is_infinite() && r.pi_ratio > 0.0, "πN>0, πS=0 must be +∞, got {}", r.pi_ratio);

        // Nothing segregating anywhere: genuinely undefined → NaN (not ∞).
        let r0 = diversity_row(&[], &[], 10.0, 10.0, 4);
        assert!(r0.pi_ratio.is_nan(), "πN=πS=0 must be NaN, got {}", r0.pi_ratio);

        // Both sides have diversity: an ordinary finite ratio.
        let rf = diversity_row(&seg, &seg, 10.0, 10.0, 4);
        assert!(rf.pi_ratio.is_finite() && rf.pi_ratio > 0.0);
    }
}

