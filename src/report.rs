//! Self-contained interactive HTML report for the VCF pN/pS analysis.
//!
//! Emits a single `.html` file with the per-gene data embedded as JSON and a
//! small vanilla-JS front-end (summary cards, an interactive Manhattan plot,
//! and a sortable/filterable table). No external assets, so it works offline.

use crate::stats;
use crate::vcf_analysis::{GenePnPs, GenomeWidePnPs};
use std::fmt::Write as _;
use std::fs::File;
use std::io::Write as _;

/// Visualization data collected from a FASTA run for the interactive report.
#[derive(Default)]
pub struct FastaReportData {
    /// (genome, lineage, dN/dS) triples from `--lineage`.
    pub lineage: Option<Vec<(String, String, f64)>>,
    /// Per-group mean dN/dS with 95% CI from `--group-average`.
    pub group: Option<Vec<crate::plot::GroupPlotData>>,
}

/// Run parameters shown in the report header.
pub struct ReportMeta<'a> {
    pub n_samples: usize,
    pub genetic_code: &'a str,
    pub kappa: f64,
    pub af_weighted: bool,
    pub fdr: f64,
    pub min_snps: usize,
    pub mk: bool,
    pub mk_fixed_af: f64,
    pub gw_ci: Option<(f64, f64)>,
    /// Genomic-control inflation factor λ (NaN if < 2 genes tested).
    pub lambda: f64,
    /// Whether the genomic-control correction was applied (p_gc/q_gc populated).
    pub genomic_control: bool,
    /// Whether repetitive genes were excluded from the pooled estimate + test.
    pub exclude_repetitive: bool,
    /// Provenance: tool version, the invoking command line, and input file paths.
    pub version: &'a str,
    pub command: &'a str,
    pub vcf_file: &'a str,
    pub ref_file: &'a str,
    pub gff_file: &'a str,
    /// Gene counts BEFORE the --min-snps filter (so "Genes analyzed"/"With SNPs"
    /// match the CLI summary and the all-genes pooled estimate, not the filtered slice).
    pub total_genes: usize,
    pub genes_with_snps: usize,
}

/// JSON-escape a string for embedding in the report. `<` and `>` are escaped as
/// `\uXXXX` too: the JSON lives inside an inline `<script>`, and an unescaped
/// `</script>` (e.g. in a gene name) would otherwise close the element early and
/// break the whole page.
fn esc(s: &str) -> String {
    let mut o = String::with_capacity(s.len() + 2);
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            '\r' => o.push_str("\\r"),
            '\t' => o.push_str("\\t"),
            '<' => o.push_str("\\u003c"),
            '>' => o.push_str("\\u003e"),
            c if (c as u32) < 0x20 => {
                let _ = write!(o, "\\u{:04x}", c as u32);
            }
            c => o.push(c),
        }
    }
    o
}

/// Format a float as a JSON literal (`null` for non-finite).
fn num(v: f64) -> String {
    if v.is_finite() {
        format!("{}", v)
    } else {
        "null".to_string()
    }
}

/// The original eskaks logo, embedded so the self-contained report needs no
/// external assets.
const LOGO_SVG: &[u8] = include_bytes!("../img/esKaKs.svg");

/// Minimal, dependency-free base64 encoder (standard alphabet, padded).
fn base64(data: &[u8]) -> String {
    const T: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(data.len().div_ceil(3) * 4);
    for chunk in data.chunks(3) {
        let b0 = chunk[0] as u32;
        let b1 = *chunk.get(1).unwrap_or(&0) as u32;
        let b2 = *chunk.get(2).unwrap_or(&0) as u32;
        let n = (b0 << 16) | (b1 << 8) | b2;
        out.push(T[(n >> 18 & 63) as usize] as char);
        out.push(T[(n >> 12 & 63) as usize] as char);
        out.push(if chunk.len() > 1 { T[(n >> 6 & 63) as usize] as char } else { '=' });
        out.push(if chunk.len() > 2 { T[(n & 63) as usize] as char } else { '=' });
    }
    out
}

/// A `data:` URI for the eskaks logo, safe to drop into an `<img src>` (the SVG's
/// own `<style>`/ids stay isolated inside the image document).
fn logo_data_uri() -> String {
    format!("data:image/svg+xml;base64,{}", base64(LOGO_SVG))
}

/// Write the interactive HTML report; returns the output path.
#[allow(clippy::too_many_arguments)]
pub fn write_html_report(
    results: &[GenePnPs],
    gw: Option<&GenomeWidePnPs>,
    core_gw: Option<&GenomeWidePnPs>,
    rep_gw: Option<&GenomeWidePnPs>,
    meta: &ReportMeta,
    divergence: Option<&std::collections::HashMap<String, f64>>,
    prefix: &str,
) -> anyhow::Result<String> {
    let output_path = format!("{}_report.html", prefix);

    // Use the PRE-filter counts (from meta) so "Genes analyzed"/"With SNPs" match the
    // CLI summary and the all-genes pooled estimate rather than the --min-snps subset.
    let total_genes = meta.total_genes;
    let genes_with_snps = meta.genes_with_snps;
    // "Tested" = the multiple-testing family (finite q); significance uses the
    // GC-corrected q under --genomic-control, matching the report's own panels.
    let n_tested = results.iter().filter(|r| r.q_value.is_finite()).count();
    let n_sig = results
        .iter()
        .filter(|r| {
            let q = if meta.genomic_control { r.q_gc } else { r.q_value };
            q.is_finite() && q < meta.fdr
        })
        .count();
    // Genome-wide SFS: sum the per-gene AF-binned counts.
    let mut sfs_n = [0u64; crate::vcf_analysis::SFS_NBINS];
    let mut sfs_s = [0u64; crate::vcf_analysis::SFS_NBINS];
    for r in results {
        for b in 0..crate::vcf_analysis::SFS_NBINS {
            sfs_n[b] += r.sfs_nonsyn[b] as u64;
            sfs_s[b] += r.sfs_syn[b] as u64;
        }
    }
    let arr = |a: &[u64]| {
        a.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(",")
    };
    let pooled_json = |g: Option<&GenomeWidePnPs>| match g {
        Some(x) => format!(
            "{{\"ratio\":{},\"pn\":{},\"ps\":{},\"nonsyn\":{},\"syn\":{}}}",
            num(x.pn_ps), num(x.pn), num(x.ps), num(x.nonsyn_snps), num(x.syn_snps)
        ),
        None => "null".to_string(),
    };

    // ── Embedded data ──────────────────────────────────────────────────────
    let mut data = String::with_capacity(4096 + results.len() * 200);
    data.push_str("{\n");

    // meta
    let ci = meta
        .gw_ci
        .map(|(lo, hi)| format!("[{}, {}]", num(lo), num(hi)))
        .unwrap_or_else(|| "null".to_string());
    let ratio_name = if meta.af_weighted { "πN/πS" } else { "pN/pS" };
    let sfs_edges = crate::vcf_analysis::SFS_EDGES
        .iter()
        .map(|e| num(*e))
        .collect::<Vec<_>>()
        .join(",");
    let _ = writeln!(
        data,
        "\"meta\":{{\"samples\":{},\"code\":\"{}\",\"kappa\":{},\"afWeighted\":{},\"fdr\":{},\"minSnps\":{},\"mk\":{},\"mkFixedAf\":{},\"ratioName\":\"{}\",\"lambda\":{},\"genomicControl\":{},\"excludeRepetitive\":{},\"version\":\"{}\",\"command\":\"{}\",\"vcfFile\":\"{}\",\"refFile\":\"{}\",\"gffFile\":\"{}\",\"sfsEdges\":[{}]}},",
        meta.n_samples, esc(meta.genetic_code), num(meta.kappa), meta.af_weighted,
        num(meta.fdr), meta.min_snps, meta.mk, num(meta.mk_fixed_af), ratio_name,
        num(meta.lambda), meta.genomic_control, meta.exclude_repetitive,
        esc(meta.version), esc(meta.command), esc(meta.vcf_file), esc(meta.ref_file),
        esc(meta.gff_file), sfs_edges
    );

    // summary
    let (gw_pn, gw_ps, gw_ratio, gw_label) = match gw {
        Some(g) => (
            num(g.pn),
            num(g.ps),
            num(g.pn_ps),
            crate::vcf_analysis::selection_label(g.pn_ps).to_string(),
        ),
        None => ("null".into(), "null".into(), "null".into(), "no data".into()),
    };
    let (gw_nsites, gw_ssites) = match gw {
        Some(g) => (num(g.n_sites), num(g.s_sites)),
        None => ("null".into(), "null".into()),
    };
    let _ = writeln!(
        data,
        "\"summary\":{{\"totalGenes\":{},\"genesWithSnps\":{},\"tested\":{},\"significant\":{},\"gwPn\":{},\"gwPs\":{},\"gwRatio\":{},\"gwLabel\":\"{}\",\"gwCi\":{},\"gwNSites\":{},\"gwSSites\":{},\"coreGw\":{},\"repGw\":{},\"sfsNonsyn\":[{}],\"sfsSyn\":[{}]}},",
        total_genes, genes_with_snps, n_tested, n_sig, gw_pn, gw_ps, gw_ratio, esc(&gw_label), ci,
        gw_nsites, gw_ssites, pooled_json(core_gw), pooled_json(rep_gw), arr(&sfs_n), arr(&sfs_s)
    );

    // genes
    data.push_str("\"genes\":[\n");
    for (i, r) in results.iter().enumerate() {
        let sites = r.n_sites + r.s_sites;
        let exp_n = if sites > 0.0 { r.n_sites / sites } else { f64::NAN };
        let comma = if i + 1 < results.len() { "," } else { "" };
        let _ = write!(
            data,
            "{{\"name\":\"{}\",\"chrom\":\"{}\",\"start\":{},\"end\":{},\"strand\":\"{}\",\"length_bp\":{},\"nSites\":{},\"sSites\":{},\"expN\":{},\"pn\":{},\"ps\":{},\"ratio\":{},\"ratioLo\":{},\"ratioHi\":{},\"nonsyn\":{},\"syn\":{},\"total\":{},\"p\":{},\"nlp\":{},\"q\":{},\"bonf\":{},\"pGc\":{},\"qGc\":{},\"rep\":{}",
            esc(&r.name), esc(&r.chrom), r.genome_start, r.genome_end, r.strand, r.length_bp,
            num(r.n_sites), num(r.s_sites), num(exp_n), num(r.pn), num(r.ps), num(r.pn_ps),
            num(r.pn_ps_lo), num(r.pn_ps_hi),
            num(r.nonsyn_snps), num(r.syn_snps), num(r.total_snps),
            num(r.p_value), num(r.neglog10p), num(r.q_value), num(r.p_bonferroni),
            num(r.p_gc), num(r.q_gc), r.repetitive
        );
        if meta.mk {
            let (dn, ds, pn, ps) =
                (r.mk_dn as f64, r.mk_ds as f64, r.mk_pn as f64, r.mk_ps as f64);
            let ni = if ps > 0.0 && dn > 0.0 { (pn * ds) / (ps * dn) } else { f64::NAN };
            let alpha = if dn > 0.0 && ps > 0.0 { 1.0 - (ds * pn) / (dn * ps) } else { f64::NAN };
            let fp = stats::fisher_exact_two_sided(
                r.mk_dn as u64, r.mk_ds as u64, r.mk_pn as u64, r.mk_ps as u64,
            );
            let _ = write!(
                data,
                ",\"dn\":{},\"ds\":{},\"pnMk\":{},\"psMk\":{},\"ni\":{},\"alpha\":{},\"fisherP\":{}",
                r.mk_dn, r.mk_ds, r.mk_pn, r.mk_ps, num(ni), num(alpha), num(fp)
            );
        }
        if let Some(div) = divergence {
            let dv = div.get(&r.name).copied().unwrap_or(f64::NAN);
            let _ = write!(data, ",\"div\":{}", num(dv));
        }
        let _ = writeln!(data, "}}{}", comma);
    }
    data.push_str("]\n}");

    // ── Assemble the HTML ──────────────────────────────────────────────────
    let mut html = String::with_capacity(HEAD.len() + BODY.len() + SCRIPT.len() + data.len() + 256);
    html.push_str(HEAD);
    html.push_str(&BODY.replace("__ESKAKS_LOGO__", &logo_data_uri()));
    html.push_str("<script>\nconst DATA = ");
    html.push_str(&data);
    html.push_str(";\n");
    html.push_str(SCRIPT);
    html.push_str("</script>\n</body>\n</html>\n");

    let mut file = File::create(&output_path)?;
    file.write_all(html.as_bytes())?;
    Ok(output_path)
}

const HEAD: &str = include_str!("report/head.html");

const BODY: &str = include_str!("report/body.html");

const SCRIPT: &str = include_str!("report/script.js");

/// Write the interactive HTML report for a FASTA (dN/dS) run: summary cards,
/// a lineage strip-scatter (points per genome + per-lineage mean), a group
/// mean±CI scatter, and the pairwise dN/dS distribution — whichever apply.
pub fn write_fasta_report(
    prefix: &str,
    model: &str,
    summary: Option<&crate::stats::SummaryStats>,
    lineage: Option<&[(String, String, f64)]>,
    group: Option<&[crate::plot::GroupPlotData]>,
    dn_ds: Option<&[(f64, f64)]>,
    window: Option<&[(usize, f64)]>,
) -> anyhow::Result<String> {
    use std::sync::atomic::Ordering;
    let output_path = format!("{}_report.html", prefix);

    // Summary values + histogram.
    let (total, valid, pooled, mean_dn, mean_ds, hist) = match summary {
        Some(s) => {
            let total = s.total_count.load(Ordering::Relaxed);
            let f = s.floats.lock().expect("summary mutex");
            let valid = f.valid_count;
            let pooled = if f.sum_ds > 0.0 { f.sum_dn / f.sum_ds } else { f64::NAN };
            let mean_dn = if valid > 0 { f.sum_dn / valid as f64 } else { f64::NAN };
            let mean_ds = if valid > 0 { f.sum_ds / valid as f64 } else { f64::NAN };
            drop(f);
            let hist = if total > 0 { Some(s.get_histogram()) } else { None };
            (total, valid, pooled, mean_dn, mean_ds, hist)
        }
        None => (0, 0, f64::NAN, f64::NAN, f64::NAN, None),
    };

    let mut data = String::with_capacity(4096);
    data.push_str("{\n");
    let _ = writeln!(
        data,
        "\"meta\":{{\"model\":\"{}\",\"totalPairs\":{},\"validPairs\":{},\"pooled\":{},\"meanDn\":{},\"meanDs\":{}}},",
        esc(model), total, valid, num(pooled), num(mean_dn), num(mean_ds)
    );

    data.push_str("\"lineage\":");
    match lineage {
        Some(lin) if !lin.is_empty() => {
            data.push('[');
            for (i, (g, l, r)) in lin.iter().enumerate() {
                let c = if i + 1 < lin.len() { "," } else { "" };
                let _ = write!(
                    data,
                    "{{\"genome\":\"{}\",\"lineage\":\"{}\",\"ratio\":{}}}{}",
                    esc(g), esc(l), num(*r), c
                );
            }
            data.push_str("],\n");
        }
        _ => data.push_str("null,\n"),
    }

    data.push_str("\"group\":");
    match group {
        Some(gr) if !gr.is_empty() => {
            data.push('[');
            for (i, g) in gr.iter().enumerate() {
                let c = if i + 1 < gr.len() { "," } else { "" };
                let _ = write!(
                    data,
                    "{{\"label\":\"{}\",\"mean\":{},\"ciLow\":{},\"ciHigh\":{}}}{}",
                    esc(&g.label), num(g.mean), num(g.ci_low), num(g.ci_high), c
                );
            }
            data.push_str("],\n");
        }
        _ => data.push_str("null,\n"),
    }

    data.push_str("\"hist\":");
    match hist {
        Some(h) if !h.is_empty() => {
            data.push('[');
            for (i, (label, count)) in h.iter().enumerate() {
                let c = if i + 1 < h.len() { "," } else { "" };
                let _ = write!(data, "{{\"label\":\"{}\",\"count\":{}}}{}", esc(label), count, c);
            }
            data.push_str("],\n");
        }
        _ => data.push_str("null,\n"),
    }

    // dN vs dS scatter (one point per pair) — a compact [dn, ds] array.
    data.push_str("\"dnds\":");
    match dn_ds {
        Some(pairs) if !pairs.is_empty() => {
            data.push('[');
            for (i, (dn, ds)) in pairs.iter().enumerate() {
                let c = if i + 1 < pairs.len() { "," } else { "" };
                let _ = write!(data, "[{},{}]{}", num(*dn), num(*ds), c);
            }
            data.push_str("],\n");
        }
        _ => data.push_str("null,\n"),
    }

    // Sliding-window dN/dS along the alignment — a positional "Manhattan".
    data.push_str("\"window\":");
    match window {
        Some(w) if !w.is_empty() => {
            data.push('[');
            for (i, (pos, r)) in w.iter().enumerate() {
                let c = if i + 1 < w.len() { "," } else { "" };
                let _ = write!(data, "[{},{}]{}", pos, num(*r), c);
            }
            data.push_str("]\n");
        }
        _ => data.push_str("null\n"),
    }
    data.push('}');

    let mut html = String::with_capacity(HEAD.len() + FASTA_BODY.len() + FASTA_SCRIPT.len() + data.len());
    html.push_str(HEAD);
    html.push_str(&FASTA_BODY.replace("__ESKAKS_LOGO__", &logo_data_uri()));
    html.push_str("<script>\nconst DATA = ");
    html.push_str(&data);
    html.push_str(";\n");
    html.push_str(FASTA_SCRIPT);
    html.push_str("</script>\n</body>\n</html>\n");

    let mut file = File::create(&output_path)?;
    file.write_all(html.as_bytes())?;
    Ok(output_path)
}

const FASTA_BODY: &str = include_str!("report/fasta_body.html");

const FASTA_SCRIPT: &str = include_str!("report/fasta_script.js");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base64_matches_rfc4648_vectors() {
        assert_eq!(base64(b""), "");
        assert_eq!(base64(b"f"), "Zg==");
        assert_eq!(base64(b"fo"), "Zm8=");
        assert_eq!(base64(b"foo"), "Zm9v");
        assert_eq!(base64(b"foob"), "Zm9vYg==");
        assert_eq!(base64(b"fooba"), "Zm9vYmE=");
        assert_eq!(base64(b"foobar"), "Zm9vYmFy");
    }

    #[test]
    fn logo_data_uri_is_a_nonempty_svg_data_uri() {
        let uri = logo_data_uri();
        assert!(uri.starts_with("data:image/svg+xml;base64,"));
        assert!(uri.len() > 1000, "logo data URI unexpectedly small");
    }
}
