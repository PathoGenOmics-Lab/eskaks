//! Per-gene pN/pS and McDonald-Kreitman table writers.

use super::*;
// Gene/chrom names come from GFF3 attributes (percent-decoded), so they can carry
// quotes, delimiters, or control chars; escape/quote them like every other writer.
use crate::textfmt::{delim_field, json_escape};

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
                    delim_field(&r.name, sep), r.length_bp, r.n_sites, r.s_sites,
                    r.pn, r.ps, format_ratio(r.pn_ps),
                    r.nonsyn_snps, r.syn_snps, r.total_snps,
                    delim_field(&r.chrom, sep), r.genome_start, r.genome_end, r.strand,
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
                        delim_field(&g.name, sep), delim_field(&g.chrom, sep), v.pos, g.strand,
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

/// Recover a segregating derived-allele count from an allele frequency and the
/// sample size, clamped to a valid segregating range `1..=n-1`.
fn derived_count(af: f64, n: usize) -> usize {
    ((af * n as f64).round() as usize).clamp(1, n - 1)
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
        for v in &g.variants {
            match v.effect {
                SnpEffect::Synonymous => syn.push(derived_count(v.af, n)),
                SnpEffect::Missense => mis.push(derived_count(v.af, n)),
                _ => {} // nonsense/stop-loss excluded, as in the pN/pS counts
            }
        }
    }
    let s_seg = syn.len() + mis.len();
    let all: Vec<usize> = syn.iter().chain(mis.iter()).copied().collect();
    let pi_n = if n_sites > 0.0 { crate::stats::theta_pi(n, &mis) / n_sites } else { f64::NAN };
    let pi_s = if s_sites > 0.0 { crate::stats::theta_pi(n, &syn) / s_sites } else { f64::NAN };
    let pi_n_pi_s = if pi_s > 0.0 { pi_n / pi_s } else { f64::NAN };
    let total_sites = n_sites + s_sites;
    let theta_w_per_site = if total_sites > 0.0 {
        crate::stats::theta_watterson(n, s_seg) / total_sites
    } else {
        f64::NAN
    };
    let tajima_d = crate::stats::tajimas_d(n, s_seg, crate::stats::theta_pi(n, &all));
    Some(GenomeDiversity { n, s_seg, pi_n, pi_s, pi_n_pi_s, theta_w_per_site, tajima_d })
}

/// Write the per-gene diversity table (`<output>_diversity.<ext>`): sample size,
/// segregating coding SNPs, per-site πN and πS (nucleotide diversity, the
/// within-species analogue of pN/pS) and their ratio, per-site Watterson θ, and
/// Tajima's D — the SFS neutrality test (D < 0: excess of rare variants; D > 0:
/// intermediate-frequency excess). Requires the sample size `n ≥ 2`.
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
    let sep = format.separator();

    writeln!(
        file,
        "Gene{s}Chrom{s}N_samples{s}S_seg{s}piN{s}piS{s}piN/piS{s}Theta_W{s}Tajima_D",
        s = sep
    )?;
    for g in results {
        let (mut syn, mut mis) = (Vec::new(), Vec::new());
        for v in &g.variants {
            match v.effect {
                SnpEffect::Synonymous => syn.push(derived_count(v.af, n)),
                SnpEffect::Missense => mis.push(derived_count(v.af, n)),
                _ => {}
            }
        }
        let s_seg = syn.len() + mis.len();
        let all: Vec<usize> = syn.iter().chain(mis.iter()).copied().collect();
        let pi_n = if g.n_sites > 0.0 { crate::stats::theta_pi(n, &mis) / g.n_sites } else { f64::NAN };
        let pi_s = if g.s_sites > 0.0 { crate::stats::theta_pi(n, &syn) / g.s_sites } else { f64::NAN };
        let pi_ratio = if pi_s > 0.0 { pi_n / pi_s } else { f64::NAN };
        let total_sites = g.n_sites + g.s_sites;
        let theta_w = if total_sites > 0.0 {
            crate::stats::theta_watterson(n, s_seg) / total_sites
        } else {
            f64::NAN
        };
        let tajima_d = crate::stats::tajimas_d(n, s_seg, crate::stats::theta_pi(n, &all));
        writeln!(
            file,
            "{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}",
            delim_field(&g.name, sep), delim_field(&g.chrom, sep), n, s_seg,
            format_pval(pi_n), format_pval(pi_s), format_ratio(pi_ratio),
            format_pval(theta_w), format_ratio(tajima_d),
            s = sep
        )?;
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
                delim_field(&r.name, sep), delim_field(&r.chrom, sep), r.genome_start, r.genome_end, r.strand,
                r.mk_dn, r.mk_ds, r.mk_pn, r.mk_ps,
                format_pval(ni(r)), format_pval(alpha(r)),
                format_pval(pvals[i]), format_pval(qvals[i]),
                s = sep
            )?;
        }
    }

    Ok(output_path)
}

