//! Per-gene pN/pS and McDonald-Kreitman table writers.

use super::*;

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
                    "  {{\"gene\":\"{}\",\"chrom\":\"{}\",\"start\":{},\"end\":{},\"strand\":\"{}\",\"length_bp\":{},\"N_sites\":{:.4},\"S_sites\":{:.4},\"exp_N_frac\":{},\"pN\":{},\"pS\":{},\"pN_pS\":{},\"nonsyn_snps\":{:.4},\"syn_snps\":{:.4},\"total_snps\":{:.4},\"p_value\":{},\"q_value_bh\":{},\"p_bonferroni\":{}}}{}",
                    r.name, r.chrom, r.genome_start, r.genome_end, r.strand, r.length_bp,
                    r.n_sites, r.s_sites, format_json_num(exp_n_frac(r)),
                    format_json_f64(r.pn), format_json_f64(r.ps), format_json_f64(r.pn_ps),
                    r.nonsyn_snps, r.syn_snps, r.total_snps,
                    format_json_num(r.p_value), format_json_num(r.q_value), format_json_num(r.p_bonferroni),
                    comma
                )?;
            }
            writeln!(file, "]")?;
        }
        _ => {
            let sep = format.separator();
            let mut file = BufWriter::new(File::create(&output_path)?);
            writeln!(
                file,
                "Gene{s}Length_bp{s}N_sites{s}S_sites{s}pN{s}pS{s}pN/pS{s}Nonsyn_SNPs{s}Syn_SNPs{s}Total_SNPs{s}Chrom{s}Start{s}End{s}Strand{s}Exp_N_frac{s}P_value{s}Q_value_BH{s}P_Bonferroni",
                s = sep
            )?;
            for r in results {
                writeln!(
                    file,
                    "{}{s}{}{s}{:.4}{s}{:.4}{s}{:.6}{s}{:.6}{s}{}{s}{:.4}{s}{:.4}{s}{:.4}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}",
                    r.name, r.length_bp, r.n_sites, r.s_sites,
                    r.pn, r.ps, format_ratio(r.pn_ps),
                    r.nonsyn_snps, r.syn_snps, r.total_snps,
                    r.chrom, r.genome_start, r.genome_end, r.strand,
                    format_pval(exp_n_frac(r)), format_pval(r.p_value),
                    format_pval(r.q_value), format_pval(r.p_bonferroni),
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
                r.name, r.chrom, r.genome_start, r.genome_end, r.strand,
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
                r.name, r.chrom, r.genome_start, r.genome_end, r.strand,
                r.mk_dn, r.mk_ds, r.mk_pn, r.mk_ps,
                format_pval(ni(r)), format_pval(alpha(r)),
                format_pval(pvals[i]), format_pval(qvals[i]),
                s = sep
            )?;
        }
    }

    Ok(output_path)
}

