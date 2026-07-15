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

    let total_genes = results.len();
    let genes_with_snps = results.iter().filter(|r| r.total_snps > 0.0).count();
    let n_tested = results.iter().filter(|r| r.p_value.is_finite()).count();
    let n_sig = results
        .iter()
        .filter(|r| r.q_value.is_finite() && r.q_value < meta.fdr)
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
    html.push_str(BODY);
    html.push_str("<script>\nconst DATA = ");
    html.push_str(&data);
    html.push_str(";\n");
    html.push_str(SCRIPT);
    html.push_str("</script>\n</body>\n</html>\n");

    let mut file = File::create(&output_path)?;
    file.write_all(html.as_bytes())?;
    Ok(output_path)
}

const HEAD: &str = r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>eskaks report</title>
<style>
/* Palette: PathoGenOmics-Lab mycolorsTB (pathogenomics brand + M. tuberculosis lineage colors). */
:root{color-scheme:light;--bg:#f9f9f7;--surface:#fcfcfb;--fg:#0b0b0b;--muted:#7d746c;
--card:#fcfcfb;--border:#e3ddd6;--grid:#e3ddd6;--axis:#c0b3a7;
--accent:#305595;--pos:#c01718;--sig:#c01718;--ns:#c0b3a7;--line:#3c5824;--sel:#d1ae00;
--s1:#ff3091;--s2:#001aff;--s3:#8a0bd2;--s4:#ff0000;--s5:#995200;--s6:#1eb040;--s7:#ff9d00;--s8:#73c2ff;}
@media (prefers-color-scheme:dark){:root:where(:not([data-theme="light"])){color-scheme:dark;
--bg:#0d0d0d;--surface:#1a1a19;--fg:#ffffff;--muted:#a89f96;--card:#1d2026;--border:#2c2c2a;--grid:#2c2c2a;--axis:#3a3a37;
--accent:#9ec4e8;--pos:#e35b5c;--sig:#e35b5c;--ns:#7d746c;--line:#7fa860;--sel:#d1ae00;}}
:root[data-theme="dark"]{color-scheme:dark;--bg:#0d0d0d;--surface:#1a1a19;--fg:#ffffff;--muted:#a89f96;
--card:#1d2026;--border:#2c2c2a;--grid:#2c2c2a;--axis:#3a3a37;
--accent:#9ec4e8;--pos:#e35b5c;--sig:#e35b5c;--ns:#7d746c;--line:#7fa860;--sel:#d1ae00;}
/* Colour-blind–safe (Okabe-Ito) palette, opt-in via the CVD toggle. Validated with
   the dataviz palette checker (adjacent-pair CVD ΔE ≥ 17); paired with shape coding. */
:root[data-cvd]{--accent:#0072B2;--pos:#D55E00;--sig:#D55E00;--ns:#999999;--line:#117733;--sel:#E69F00;
--s1:#0072B2;--s2:#D55E00;--s3:#CC79A7;--s4:#009E73;--s5:#E69F00;--s6:#56B4E9;--s7:#F0E442;--s8:#117733;}
:root[data-cvd][data-theme="dark"]{--ns:#8a8f98;}
*{box-sizing:border-box}
body{margin:0;font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;
background:var(--bg);color:var(--fg);line-height:1.5}
.wrap{max-width:1280px;margin:0 auto;padding:24px}
.layout{display:flex;gap:20px;align-items:flex-start}
.main{flex:1;min-width:0}
.toc{position:sticky;top:16px;flex:0 0 196px;max-height:calc(100vh - 32px);overflow:auto;
border:1px solid var(--border);border-radius:10px;background:var(--surface);padding:10px 6px}
.toc .toc-h{font-size:.7rem;text-transform:uppercase;letter-spacing:.05em;color:var(--muted);padding:2px 8px 6px}
.toc a{display:block;padding:4px 8px;font-size:.8rem;line-height:1.3;color:var(--muted);text-decoration:none;
border-left:2px solid transparent;border-radius:0 5px 5px 0;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.toc a:hover{color:var(--fg)}
.toc a.active{color:var(--accent);border-left-color:var(--accent);background:var(--card);font-weight:600}
@media (max-width:900px){.toc{display:none}.layout{display:block}}
@media print{.toc{display:none}.layout{display:block}}
h1{font-size:1.5rem;margin:0 0 4px}
.sub{color:var(--muted);font-size:.85rem;margin-bottom:14px}
.topbar{display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-bottom:18px}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin-bottom:20px}
.card{background:var(--card);border:1px solid var(--border);border-radius:10px;padding:14px}
.card .k{color:var(--muted);font-size:.72rem;text-transform:uppercase;letter-spacing:.04em}
.card .v{font-size:1.3rem;font-weight:650;margin-top:4px}
.card .v.small{font-size:.98rem}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:20px;margin-bottom:24px}
section{margin-bottom:8px}
.panel{background:var(--surface);border:1px solid var(--border);border-radius:12px;padding:14px}
.panel.wide{grid-column:1 / -1}
h2{font-size:1rem;margin:0 0 8px;padding-bottom:6px;border-bottom:1px solid var(--border);
display:flex;align-items:center;gap:8px}
h2 .htitle{display:inline-flex;align-items:center;gap:8px;min-width:0;flex:1 1 auto}
.toolbar{display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin:8px 0}
button.tog,button.btn{background:var(--card);border:1px solid var(--border);color:var(--fg);
padding:6px 11px;border-radius:8px;cursor:pointer;font-size:.83rem}
button.tog.on{background:var(--accent);color:#fff;border-color:var(--accent)}
button.btn:hover,button.tog:hover{border-color:var(--accent)}
input[type=search]{flex:1;min-width:180px;padding:7px 10px;border:1px solid var(--border);
border-radius:8px;background:var(--surface);color:var(--fg);font-size:.9rem}
.legend{display:flex;gap:12px;flex-wrap:wrap;font-size:.75rem;color:var(--muted);margin-top:4px}
.legend span{display:inline-flex;align-items:center;gap:5px}
.legend i{width:11px;height:11px;border-radius:50%;display:inline-block}
svg{max-width:100%;height:auto;display:block}
.tip{position:fixed;pointer-events:none;background:var(--fg);color:var(--bg);
padding:6px 9px;border-radius:6px;font-size:.78rem;opacity:0;transition:opacity .08s;z-index:20;white-space:nowrap;box-shadow:0 2px 8px rgba(0,0,0,.25)}
.tablewrap{overflow-x:auto;border:1px solid var(--border);border-radius:10px;max-height:560px}
table{border-collapse:collapse;width:100%;font-size:.82rem;font-variant-numeric:tabular-nums}
th,td{padding:6px 9px;text-align:right;white-space:nowrap}
th:first-child,td:first-child{text-align:left;position:sticky;left:0;background:var(--surface)}
thead th{position:sticky;top:0;background:var(--card);cursor:pointer;user-select:none;border-bottom:1px solid var(--border)}
thead th:hover{color:var(--accent)}
thead th.sorted::after{content:" ▲";font-size:.7em}
thead th.sorted.desc::after{content:" ▼"}
tbody tr.zebra{background:var(--card)}
tbody tr.sig td{font-weight:600}
tbody tr.sel{outline:2px solid var(--sel);outline-offset:-2px}
tbody tr{cursor:pointer}
.badge{display:inline-block;font-size:.62rem;padding:1px 5px;border-radius:5px;background:var(--pos);color:#fff;margin-left:4px}
.mbar{display:inline-flex;height:9px;width:70px;border-radius:3px;overflow:hidden;vertical-align:middle;border:1px solid var(--border)}
.count{color:var(--muted);font-size:.8rem;margin-left:auto}
details.methods{margin:8px 0 20px;border:1px solid var(--border);border-radius:10px;background:var(--card)}
details.methods summary{cursor:pointer;padding:10px 14px;font-weight:600}
details.methods .body{padding:0 14px 12px;font-size:.82rem;color:var(--muted);display:flex;gap:8px;flex-wrap:wrap}
.prov{flex-basis:100%;margin-top:8px;border-top:1px dashed var(--border);padding-top:8px;line-height:1.7}
.prov code{background:var(--surface);border:1px solid var(--border);border-radius:4px;padding:1px 5px;font-size:.92em;word-break:break-all}
.chip{background:var(--surface);border:1px solid var(--border);border-radius:20px;padding:3px 10px}
.census{display:flex;gap:8px;flex-wrap:wrap}
.chip.regime{cursor:pointer;color:var(--fg);display:inline-flex;align-items:center;gap:6px;font-size:.82rem}
.chip.regime i{width:11px;height:11px;border-radius:50%;display:inline-block}
.chip.regime.on{border-color:var(--accent);background:var(--card)}
.chip.regime b{font-variant-numeric:tabular-nums}
.candlist{margin-top:8px;font-size:.82rem;display:flex;gap:6px;flex-wrap:wrap;align-items:center}
.cand{cursor:pointer;background:var(--card);border:1px solid var(--border);border-radius:6px;padding:2px 7px}
.cand:hover{border-color:var(--sel)}
.cand small{color:var(--muted)}
footer{color:var(--muted);font-size:.75rem;margin-top:24px;text-align:center}
.muted{color:var(--muted);font-size:.85rem}
.visually-hidden{position:absolute;width:1px;height:1px;overflow:hidden;clip:rect(0 0 0 0);white-space:nowrap}
circle.mark,rect.mark,path.mark{cursor:pointer}
/* ── "i" interpretation help ────────────────────────────────────── */
.info{width:19px;height:19px;min-width:19px;border-radius:50%;border:1px solid var(--border);
background:var(--card);color:var(--muted);font:italic 700 .72rem/1 Georgia,serif;cursor:pointer;
display:inline-flex;align-items:center;justify-content:center;padding:0;flex:none}
.info:hover,.info.on{border-color:var(--accent);color:var(--accent)}
.card .info{width:16px;height:16px;min-width:16px;vertical-align:middle;margin-left:5px}
.info-pop{position:fixed;z-index:40;max-width:340px;background:var(--card);border:1px solid var(--border);
border-radius:10px;padding:12px 14px;box-shadow:0 8px 28px rgba(0,0,0,.28);font-size:.8rem;color:var(--fg);line-height:1.45}
.info-pop h4{margin:0 0 5px;font-size:.85rem;color:var(--accent)}
.info-pop .lead{margin:0 0 8px}
.info-pop .ax{color:var(--muted);font-size:.75rem;margin:2px 0}
.info-pop .ax b{color:var(--fg);font-style:normal}
.info-pop ul{margin:7px 0 0;padding-left:16px}
.info-pop li{margin:3px 0}
.info-pop .warn{margin-top:9px;color:var(--muted);font-size:.75rem;border-top:1px dashed var(--border);padding-top:7px}
.info-pop .warn b{color:var(--sel);font-style:normal}
.guide{margin:8px 0 20px;border:1px solid var(--border);border-radius:10px;background:var(--card)}
.guide summary{cursor:pointer;padding:10px 14px;font-weight:600}
.guide .gbody{padding:0 14px 14px;font-size:.83rem;line-height:1.5}
.guide .gloss{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:10px 20px;margin-top:6px}
.guide .gloss div{color:var(--muted)}
.guide .gloss b{color:var(--fg)}
.verdict{display:flex;gap:14px;align-items:center;flex-wrap:wrap;padding:12px 16px;border-radius:10px;
border:1px solid var(--border);background:var(--card);margin-bottom:18px}
.verdict .vlab{font-size:1.05rem;font-weight:650}
.verdict .vsub{color:var(--muted);font-size:.83rem}
.verdict .dot{width:12px;height:12px;border-radius:50%;flex:none}
@media print{.topbar,.info,.guide summary::marker{display:none!important}.panel,.verdict{break-inside:avoid}
.grid{display:block}.panel{margin-bottom:16px}.guide[open] summary{display:none}}
</style>
</head>
<body>
"#;

const BODY: &str = r#"<div class="wrap">
<h1>eskaks — pN/pS report</h1>
<div class="sub" id="meta"></div>

<div class="topbar">
  <input type="search" id="search" placeholder="Search gene / chromosome…">
  <button class="tog" id="strTog" title="Multiple-testing stringency">FDR (BH)</button>
  <button class="btn" id="expCsv">⤓ CSV</button>
  <button class="btn" id="expJson">⤓ JSON</button>
  <button class="btn" id="printBtn" title="Print or save as PDF">🖶 Print</button>
  <button class="btn" id="cvdTog" title="Colour-blind–safe palette + shapes" aria-pressed="false">👁 CVD</button>
  <button class="btn" id="themeTog" title="Toggle light/dark">◑ Theme</button>
</div>

<div class="layout">
<aside class="toc" id="toc" aria-label="Table of contents"></aside>
<main class="main">
<div class="verdict" id="verdict"></div>

<div class="cards" id="cards"></div>

<details class="guide"><summary>▸ How to read this report</summary><div class="gbody" id="guide"></div></details>
<details class="methods"><summary>Methods &amp; parameters</summary><div class="body" id="methods"></div></details>

<div class="grid" id="panels"></div>

<section id="sec-table">
<h2>Per-gene table <span class="count" id="tableCount" aria-live="polite"></span></h2>
<div class="tablewrap" id="tblwrap"><table id="tbl"><caption class="visually-hidden">Per-gene pN/pS results; click a row to highlight that gene in every panel. Use up and down arrows to step through genes.</caption><thead></thead><tbody></tbody></table></div>
</section>
</main>
</div>

<footer>Generated by eskaks · self-contained, no external assets · click a point or row to highlight a gene everywhere</footer>
</div>
<div class="tip" id="tip"></div>
<div class="info-pop" id="infopop" hidden></div>
"#;

const SCRIPT: &str = r##"
const $ = s => document.querySelector(s);
const genes = DATA.genes, S = DATA.summary, M = DATA.meta;
const fmt = (v,d=4) => v==null||!isFinite(v) ? "NA" : Number(v).toFixed(d);
const fmtP = v => v==null||!isFinite(v) ? "NA" : (v!==0 && Math.abs(v)<1e-3 ? v.toExponential(2) : v.toFixed(4));
const tip = $("#tip");
genes.forEach((g,i)=>{ g._i=i;
  // Power-aware effect size: binomial z of observed nonsyn SNPs vs the neutral expectation total*expN.
  const e=(g.expN!=null&&isFinite(g.expN))?g.expN:null;
  g.z=(e!=null&&e>0&&e<1&&g.total>0)?(g.nonsyn-g.total*e)/Math.sqrt(g.total*e*(1-e)):null;
  // Direction of Selection (MK): stays finite at low counts, unlike alpha.
  if(g.dn!=null){ const dd=g.dn+g.ds, pp=g.pnMk+g.psMk; g.dos=(dd>0&&pp>0)?(g.dn/dd)-(g.pnMk/pp):null; }
});

// ── State ─────────────────────────────────────────────────────────────
let metric = "neglogp";      // Manhattan y-axis
let stringency = "q";        // q (BH) | bonf (Bonferroni)
let selected = null;         // selected gene index
let cvd = false;             // colour-blind–safe palette + shape coding
const RE_QUAR = /(PE_PGRS|PPE|^PE|PE\d|PPE\d|PGRS|maturase|transpos|IS6110)/i;
// Significance value at the current stringency. When genomic control was applied,
// the BH branch uses the GC-corrected q so the whole report reflects the correction.
const sigVal = g => stringency==="q" ? (M.genomicControl ? g.qGc : g.q) : g.bonf;
// The p-value the significance boundary is built on: GC-corrected under
// --genomic-control (BH mode), else the raw test p. Keeps the Manhattan/volcano
// y-axis and threshold line consistent with the qGc-based colouring.
const pStat = g => (M.genomicControl && stringency==="q") ? g.pGc : g.p;
const isSig  = g => { const v=sigVal(g); return v!=null && isFinite(v) && v < M.fdr; };
const quar   = g => (g.rep!=null ? g.rep : RE_QUAR.test(g.name||""));
// Is the gene diversifying (pN/pS>1)? When pN/pS is +∞ (syn=0, nonsyn>0) the ratio
// is serialized as null, so fall back to the counts — that gene is the STRONGEST
// diversifying signal and must not be mislabelled purifying.
const isUp = g => (g.ratio!=null&&isFinite(g.ratio)) ? g.ratio>1 : ((g.nonsyn||0) > (g.syn||0));
// Direction as a marker shape (used in CVD mode so colour is never the only cue).
const dirShape = g => !isSig(g) ? "dot" : (isUp(g) ? "up" : "down");
// Selection regime from polymorphism + significance. Genes that fail the test are
// "not significant" (underpowered) — NOT asserted neutral: they may still be under
// selection but lack the SNPs to reject H0.
function regime(g){ if(!isSig(g)) return "ns"; return isUp(g) ? "positive" : "purifying"; }
const REGIMES = {positive:"var(--pos)", purifying:"var(--accent)", ns:"var(--ns)"};
const RLAB = {positive:"positive", purifying:"purifying", ns:"not significant"};
let regimeFilter = null;   // regime name to filter the table by, or null

// ── "i" interpretation help ───────────────────────────────────────────
const HELP = {"census":{"title":"Selection regimes census","lead":"Tallies genes into positive, purifying, or not-significant selection based on their pN/pS ratio and neutrality-test outcome.","x":"","y":"","read":["Positive chip: significant with pN/pS&gt;1 — excess amino-acid changes; scan for drug-resistance or immune-epitope loci.","Purifying chip: significant with pN/pS&lt;1 — the genome-wide baseline for essential, conserved M. tuberculosis genes.","Not-significant chip: neutrality not rejected, usually too few SNPs; do NOT read as evidence of neutrality.","Click any chip to filter the per-gene table and inspect the genes driving that count."],"watch":["Most genes carry few SNPs in clonal, low-diversity M. tuberculosis, so 'not significant' reflects low power, not neutrality.","PE/PPE and repetitive regions give mapping artefacts and spurious SNPs; treat their positive calls sceptically. Use FDR q, not raw p."]},"manhattan":{"title":"Manhattan: per-gene selection scan","lead":"Each gene plotted by genome position against the strength of its departure from neutral pN/pS, revealing which loci are under selection and where they sit.","x":"Genome position in bp along the H37Rv reference; genes ordered left-to-right, so vertical stripes cluster co-located loci.","y":"-log10(p) of the neutrality test (higher = stronger evidence against pN/pS=1); toolbar toggle switches y to raw pN/pS (1 = neutral).","read":["Coloured points above the dashed line pass FDR/Bonferroni — red = diversifying, blue = purifying; candidates worth follow-up.","Colour shows direction of significant genes: red = pN/pS&gt;1 (diversifying), blue = &lt;1 (purifying); grey = not significant.","Large points = many SNPs, so the test has power; trust these over tiny points at the same height.","Isolated tall spikes flag single-gene signals; clustered peaks may reflect a locus or region."],"watch":["Low-count genes give unstable ratios and weak power; a small point just over the line is fragile.","M. tuberculosis clonality means genome-wide linkage — nearby peaks are often correlated, not independent hits; PE/PPE peaks are frequently mapping artefacts."]},"volcano":{"title":"Volcano: effect vs significance","lead":"Each gene plotted by selection direction and strength (x) against statistical confidence (y), so the corners flag the strongest, most reliable departures from neutrality.","x":"log2(pN/pS): negative = purifying selection, 0 = neutral (mutational expectation), positive = diversifying/positive selection; each unit is a twofold change.","y":"-log10(p) from the exact binomial neutrality test; higher = stronger evidence against pN/pS=1. Above the dashed line clears the significance threshold.","read":["Upper-left (x&lt;0, high y): genes under significant purifying selection, the expected signal for conserved, essential loci.","Upper-right (x&gt;0, high y): significant excess of nonsynonymous SNPs, candidate positive/diversifying selection, e.g. drug targets or antigens.","Bottom band (low y): non-significant regardless of x; ratio scatter here is noise, not selection.","Read direction (x) and confidence (y) jointly; large log2 alone means nothing if the point sits low."],"watch":["Low-count genes give extreme log2(pN/pS) at low y; confirm against the significance line, and prefer q/Bonferroni over raw p given genome-wide testing.","M. tuberculosis clonality links sites, so a driver SNP drags hitchhikers; PE/PPE points may be mapping artefacts, not real diversifying selection."]},"mk":{"title":"McDonald-Kreitman: α vs significance","lead":"Per-gene test contrasting nonsynonymous vs synonymous divergence and polymorphism to detect adaptive protein evolution.","x":"α = 1 − (Ds·Pn)/(Dn·Ps): estimated fraction of nonsynonymous divergence fixed by positive selection (dimensionless, capped near 1; can be negative).","y":"−log10(Fisher exact p) on the 2×2 Dn/Ds vs Pn/Ps table; higher = stronger evidence against neutrality.","read":["Top-right (α&gt;0, high y, red): adaptive evolution — excess fixed nonsynonymous changes; NI&lt;1. Prime positive-selection candidates.","Top-left (α&lt;0, high y, blue): excess nonsynonymous polymorphism, NI&gt;1 — segregating slightly-deleterious variants or purifying constraint.","Bottom band (low y): non-significant regardless of α; too few counts to reject neutrality — treat α as noisy.","Rank the red, high-y genes; cross-check whether they are known drug-resistance or antigen loci."],"watch":["Genes with any zero cell (Dn,Ds,Pn,Ps) give unstable α and inflated |α|; require all four counts non-trivial before believing the sign.","M. tuberculosis clonality creates genome-wide linkage, so per-gene Fisher tests are non-independent — apply FDR and distrust PE/PPE hits from repetitive-region mismapping."]},"recon":{"title":"Polymorphism vs divergence","lead":"Contrasts within-sample selection (pN/pS) against long-term divergence (dN/dS) per gene to spot genes whose selective regime has changed over time.","x":"pN/pS: nonsynonymous/synonymous polymorphism in this VCF, site-normalized; &gt;1 excess replacement, &lt;1 purifying, =1 neutral.","y":"dN/dS: nonsynonymous/synonymous substitution rate between lineages/species (supplied); &gt;1 past positive selection, &lt;1 constraint.","read":["Top-left (dN/dS&gt;1, pN/pS&lt;1): PAST-POSITIVE — historically adaptive, now purifying; enlarged clickable candidates for antigen/resistance loci.","Top-right (both&gt;1): DIVERSIFYING — sustained diversifying selection across timescales.","Bottom-right (pN/pS&gt;1, dN/dS&lt;1): RELAXED/RECENT — recent relaxation or emerging adaptation not yet fixed.","Bottom-left (both&lt;1) / on-diagonal: PURIFYING or concordant regimes; off-diagonal points flag a shifted selective pressure."],"watch":["Low SNP counts make pN/pS unstable; M. tuberculosis clonality and near-zero within-host diversity leave many genes with few/no polymorphisms.","PE/PPE and repetitive genes give unreliable dN/dS and mapping artefacts; dN/dS must share pN/pS's kappa/site model or the diagonal is misaligned."]},"funnel":{"title":"Power funnel: pN/pS vs SNP count","lead":"Each gene's pN/pS plotted against how many SNPs support it, showing which pN/pS values are statistically resolvable versus small-sample noise.","x":"log10 of total SNPs in the gene (mutational evidence; each step right is 10x more variants and a tighter estimate)","y":"Gene pN/pS ratio; 1 = neutral expectation, &lt;1 = purifying selection, &gt;1 = diversifying/positive selection","read":["Wide left mouth: few-SNP genes scatter far from 1 by chance alone; extreme pN/pS here is usually not significant.","Narrow right neck: many-SNP genes converge on the pooled dashed line, the genome's purifying baseline (typically &lt;1 in Mtb).","High-count genes sitting well above 1: strong positive-selection candidates worth the neutrality test and q-value.","Points left of the --min-snps line are underpowered and excluded from ranking; treat their pN/pS as indicative only."],"watch":["Mtb clonality and tight linkage mean SNP counts are not independent draws, so per-gene power is lower than the count implies.","PE/PPE genes inflate to the right via repetitive-region mapping/alignment artefacts, mimicking real diversifying selection."]},"obsexp":{"title":"Observed vs expected N-fraction","lead":"Plots each gene's observed nonsynonymous SNP fraction against the neutral expectation, so deviation from the diagonal reveals selection.","x":"Expected nonsynonymous fraction N/(N+S) under neutrality (0-1), set by codon composition and the transition/transversion bias kappa.","y":"Observed nonsynonymous / total SNPs in the gene (0-1); the empirical fraction of variant sites that change amino acids.","read":["On the dashed diagonal: observed matches neutral expectation, pN/pS approximately 1.","Above the line (red): excess nonsynonymous SNPs, diversifying or relaxed selection; candidate antigen/drug-target genes.","Below the line (blue): nonsynonymous deficit, purifying selection preserving protein sequence; typical for core essentials.","Grey points are non-significant; treat their apparent offset as noise, not selection."],"watch":["Low-count genes scatter widely off the diagonal by chance; a colour survives only after BH/Bonferroni, so judge significance, not raw distance.","M. tuberculosis clonality links sites, so hitchhiking can mimic per-gene selection; PE/PPE mapping artefacts and a kappa-dependent x-axis further bias placement."]},"dist":{"title":"pN/pS distribution","lead":"Genome-wide histogram of per-gene pN/pS, revealing the dominant purifying-selection peak below 1 and a right tail of candidate diversifying genes.","x":"Per-gene pN/pS: nonsynonymous-to-synonymous polymorphism, normalized by N/S mutational opportunity (1.0 = neutral)","y":"Number of genes falling in each pN/pS bin (count)","read":["Mass below 1: purifying selection removing nonsynonymous variants — the expected bacterial baseline.","Neutral bin at 1.0: relaxed constraint or too few SNPs to distinguish from neutrality.","Right tail beyond 1: candidate diversifying/positive selection — inspect these genes, but confirm with the per-gene test.","A heavy or shifted distribution overall can signal a demographic or filtering artefact, not selection."],"watch":["Low-SNP genes give unstable, extreme ratios; a bin far from 1 needs the exact test's q-value, not visual position.","M. tuberculosis clonality plus PE/PPE mapping artefacts inflate spurious nonsynonymous counts, padding the right tail."]},"qq":{"title":"P-value QQ — observed vs null","lead":"Diagnostic comparing the whole set of neutrality p-values against what pure chance (a uniform null) would produce.","x":"expected &minus;log10(p) if every gene were neutral (order statistics of a uniform distribution)","y":"observed &minus;log10(p) for the same rank; the &lambda; label is the genomic-inflation factor (median obs / median exp)","read":["Points on the diagonal (y=x): the genes behave as the neutral null predicts.","An upward-bending tail at the right: more strong signals than chance — genuine selection candidates.","&lambda; near 1 is well-calibrated; &lambda; &gt;&gt; 1 means far more low p-values than chance — widespread selection or systematic bias (mapping, reference, filtering)."],"watch":["M. tuberculosis clonality correlates genes, so an inflated tail can be linkage rather than many independent targets."]},"card_gw":{"title":"Genome-wide pN/pS (pooled)","lead":"One overall ratio pooled across all genes: &Sigma; nonsynonymous SNPs / &Sigma; N sites over &Sigma; synonymous / &Sigma; S sites.","x":"","y":"","read":["Below 1: pervasive purifying selection — the expected bacterial baseline.","The bracket is a bootstrap 95% CI from resampling genes; if it excludes 1, the departure is well supported."],"watch":["Pooling is dominated by high-SNP genes; it summarises the genome, not any single locus."]},"card_tested":{"title":"Genes tested","lead":"Genes with enough SNPs to run the neutrality test (after --min-snps), i.e. the multiple-testing family size.","x":"","y":"","read":["This is the denominator for BH-FDR and Bonferroni correction across genes."],"watch":["Many M. tuberculosis genes are untested for lack of SNPs — absence of a hit is not evidence of neutrality."]},"card_sig":{"title":"Significant genes","lead":"Genes rejecting neutrality (pN/pS&ne;1) at the current stringency — toggle FDR(BH) &harr; Bonferroni in the top bar.","x":"","y":"","read":["BH-FDR controls the expected false-discovery proportion — more permissive, good for screening.","Bonferroni controls the family-wise error — stricter, good for a confident shortlist."],"watch":["The count changes with the toggle; report which one you used."]}};
HELP.card_lambda={title:"Genomic inflation λ",lead:"How far the whole set of neutrality p-values departs from a uniform null — a diagnostic, not automatically a correction target.",x:"",y:"",
  read:["λ ≈ 1: the tested genes behave as chance predicts.","λ ≫ 1: many more low p-values than chance — this can be genuine, pervasive purifying selection OR systematic inflation (mapping, reference, filtering, structure).","--genomic-control divides each χ² by λ and re-tests (q reads 'corrected') — only appropriate when the high λ is artefactual, not real genome-wide selection."],
  watch:["In clonal M. tuberculosis, real pervasive selection routinely drives λ well above 1, so do NOT apply genomic control reflexively — reserve it for suspected systematic bias."]};
function helpHtml(h){
  return '<h4>'+h.title+'</h4><p class="lead">'+h.lead+'</p>'
    +(h.x?'<div class="ax"><b>x</b> — '+h.x+'</div>':'')
    +(h.y?'<div class="ax"><b>y</b> — '+h.y+'</div>':'')
    +((h.read&&h.read.length)?'<ul>'+h.read.map(function(s){return '<li>'+s+'</li>';}).join('')+'</ul>':'')
    +((h.watch&&h.watch.length)?'<div class="warn"><b>Watch out</b> — '+h.watch.join(' · ')+'</div>':'');
}
const infopop=$("#infopop"); let openInfo=null;
function positionPop(btn){
  const r=btn.getBoundingClientRect(), pw=infopop.offsetWidth, ph=infopop.offsetHeight, m=8;
  const left=Math.max(m, Math.min(r.left, window.innerWidth-pw-m));
  let top=r.bottom+6; if(top+ph>window.innerHeight-m) top=r.top-ph-6;   // flip above if it would overflow below
  top=Math.max(m, Math.min(top, window.innerHeight-ph-m));               // clamp fully into the viewport
  infopop.style.left=left+'px'; infopop.style.top=top+'px';
}
function showInfo(btn){ const h=HELP[btn.dataset.help]; if(!h) return;
  infopop.innerHTML=helpHtml(h); infopop.hidden=false; positionPop(btn);
  btn.classList.add('on'); openInfo=btn; }
function hideInfo(){ infopop.hidden=true; if(openInfo)openInfo.classList.remove('on'); openInfo=null; }
document.addEventListener('click', function(e){ const b=e.target.closest&&e.target.closest('.info');
  if(b){ e.stopPropagation(); if(openInfo===b){hideInfo();}else{hideInfo();showInfo(b);} return; }
  if(!(e.target.closest&&e.target.closest('#infopop'))) hideInfo(); });
document.addEventListener('keydown', function(e){ if(e.key==='Escape' && openInfo){ hideInfo(); e.stopImmediatePropagation(); } });  // 1st Esc closes the popover only
window.addEventListener('resize', hideInfo);
window.addEventListener('scroll', hideInfo, true);   // capture: also fires for the table's own scroll box

// ── Genome-wide verdict banner + glossary ─────────────────────────────
function renderVerdict(){
  const c={positive:0,purifying:0,ns:0}; genes.forEach(function(g){c[regime(g)]++;});
  const gw=S.gwRatio, lab=(S.gwLabel||'').split(' (')[0];
  // Dot colour follows the label (same thresholds as the Selection card), never independent cutoffs.
  const dotc=/purifying/i.test(lab)?'var(--accent)':(/diversif|positive/i.test(lab)?'var(--pos)':'var(--ns)');
  const ciTxt=S.gwCi?' ['+fmt(S.gwCi[0],2)+'–'+fmt(S.gwCi[1],2)+']':'';
  const strLab=stringency==='q'?'BH-FDR':'Bonferroni';
  const untested=Math.max(0,S.totalGenes-S.tested);
  $("#verdict").innerHTML=
    '<span class="dot" style="background:'+dotc+'"></span>'+
    '<div style="flex:1;min-width:0"><div class="vlab">Genome-wide '+M.ratioName+' '+fmt(gw,3)+ciTxt+' — '+lab+'</div>'+
    '<div class="vsub">'+c.positive+' diversifying · '+c.purifying+' purifying — significant ('+strLab+' &lt; '+M.fdr+') of '+S.tested+' gene(s) tested'+(untested?'; '+untested+' untested (too few SNPs)':'')+'</div></div>'+
    '<button class="info" data-help="card_gw" title="What is this?" aria-label="What is this?">i</button>';
}
function renderGuide(){
  const terms=[
    ['pN/pS','nonsynonymous vs synonymous polymorphism, each divided by its number of sites. &lt;1 purifying · ≈1 neutral · &gt;1 diversifying.'],
    ['N / S sites','the mutational opportunity for nonsynonymous / synonymous change in a gene, κ-weighted for the ts/tv bias.'],
    ['expN = N/(N+S)','fraction of random mutations expected to be nonsynonymous under neutrality — the null for the test.'],
    ['neutrality test','two-sided exact binomial of observed nonsynonymous SNPs vs expN (H0: pN/pS=1); reported as p, q(BH), p(Bonferroni).'],
    ['genome-wide pN/pS','one ratio pooled over all genes (Σ SNPs / Σ sites), with a bootstrap 95% CI.'],
    ['selection regime','significant &amp; pN/pS&gt;1 = positive · significant &amp; &lt;1 = purifying · else not significant (often underpowered).'],
    ['κ (kappa)','transition/transversion rate ratio used when counting sites; κ&gt;1 raises S (and lowers N), which raises pN/pS under a ts-skewed spectrum.'],
    ['z(N)','standardized nonsynonymous excess (obs−exp)/SD — a power-aware effect size that stays finite at low counts.'],
    ['pN/pS 95% CI','Wilson interval on the nonsynonymous SNP fraction, mapped to pN/pS; if it excludes 1 the gene departs neutrality.'],
    ['genomic inflation λ','median χ² of the tested genes ÷ its null — λ≈1 is calibrated; λ≫1 can be real pervasive selection or systematic bias.'],
    ['SFS pN/pS','pN/pS split by allele frequency; purifying selection makes it fall as frequency rises (deleterious nonsyn stay rare).'],
    ['core vs repetitive','pooled pN/pS for core genes vs PE/PPE/PGRS/IS; a big gap warns of repeat mapping artefacts.'],
  ];
  if(M.mk){ terms.push(
    ['McDonald-Kreitman','fixed (Dn,Ds) vs polymorphic (Pn,Ps) changes; NI=(Pn/Ps)/(Dn/Ds), α=1−NI, Fisher exact p.'],
    ['DoS','Direction of Selection = Dn/(Dn+Ds) − Pn/(Pn+Ps); &gt;0 adaptive, &lt;0 constrained — robust at low counts.']); }
  if(M.genomicControl){ terms.push(
    ['genomic control','the neutrality χ² was divided by λ and re-tested; the q shown is GC-corrected (p_gc/q_gc columns).']); }
  if(hasDiv){ terms.push(
    ['dN/dS (divergence)','between-lineage substitution ratio supplied via --divergence; compared against within-sample pN/pS.']); }
  $("#guide").innerHTML=
    '<p>Every gene is tested for departure from neutral evolution from its synonymous vs nonsynonymous SNPs. '+
    'Click the round <b>i</b> on any panel or card for how to read it; click a point or table row to highlight that gene across every panel; use ↑/↓ to step through genes and Esc to clear.</p>'+
    '<div class="gloss">'+terms.map(function(t){return '<div><b>'+t[0]+'</b> — '+t[1]+'</div>';}).join('')+'</div>';
}

// ── Meta + cards + methods ────────────────────────────────────────────
$("#meta").textContent =
  `${M.samples} sample(s) · genetic code ${M.code} · kappa ${fmt(M.kappa,2)}` +
  (M.afWeighted ? " · AF-weighted (πN/πS)" : "") + (M.minSnps>0 ? ` · min-snps ${M.minSnps}` : "");
const infoBtn = key => (HELP[key]?`<button class="info" data-help="${key}" title="What is this?" aria-label="What is this?">i</button>`:"");
const card=(k,v,s,help)=>`<div class="card"><div class="k">${k}${help?infoBtn(help):""}</div><div class="v${s?" small":""}">${v}</div></div>`;
const ciTxt = S.gwCi ? ` [${fmt(S.gwCi[0])}, ${fmt(S.gwCi[1])}]` : "";
const lamText = (M.lambda!=null&&isFinite(M.lambda)) ? fmt(M.lambda,2) : "NA";
function renderCards(){
  const nSig = genes.filter(isSig).length;
  const sigLbl = stringency==="q" ? (M.genomicControl?"BH·GC":"BH") : "Bonf";
  $("#cards").innerHTML =
    card("Genes analyzed", S.totalGenes) + card("With SNPs", S.genesWithSnps) +
    card(`Genome-wide ${M.ratioName}`, fmt(S.gwRatio,4)+`<span style="font-size:.6em;color:var(--muted)">${ciTxt}</span>`, true, "card_gw") +
    card("Selection", S.gwLabel.split(" (")[0], true) +
    card("Genes tested", S.tested, false, "card_tested") +
    card(`Significant (${sigLbl}<${M.fdr})`, nSig, false, "card_sig") +
    card("Inflation λ", lamText+(M.genomicControl?' <span style="font-size:.55em;color:var(--muted)">corrected</span>':''), true, "card_lambda");
}
$("#methods").innerHTML = [
  ["samples",M.samples],["genetic code",M.code],["kappa",fmt(M.kappa,2)],
  ["af-weighted",M.afWeighted],["FDR",M.fdr],["min-snps",M.minSnps],
  M.mk?["mk-fixed-af",M.mkFixedAf]:null,
  ["inflation λ",lamText],
  M.genomicControl?["genomic-control","on"]:null,
  M.excludeRepetitive?["exclude-repetitive","on"]:null,
  ["eskaks",M.version]
].filter(Boolean).map(([k,v])=>`<span class="chip">${k}: <b>${v}</b></span>`).join("")
  + `<div class="prov">Command: <code>${(M.command||"").replace(/</g,"&lt;")}</code><br>`
  + `Inputs: VCF <code>${M.vcfFile||"?"}</code> · ref <code>${M.refFile||"?"}</code> · GFF <code>${M.gffFile||"?"}</code></div>`;

// ── BH threshold p* (for the given stringency) ────────────────────────
function pThreshold(){
  const ps = genes.map(pStat).filter(v=>v!=null&&isFinite(v)).sort((a,b)=>a-b);
  const m = ps.length; if(!m) return null;
  if(stringency==="bonf") return M.fdr/m;
  let pStar=null; for(let i=0;i<m;i++){ if(ps[i] <= ((i+1)/m)*M.fdr) pStar=ps[i]; } return pStar;
}

// ── Generic scatter renderer (returns SVG string, or SVG+canvas for large N) ──
// cfg: {rows,x,y,xlabel,ylabel,xfmt,yfmt,color,size,shape,yerr,tip,refs,W,H}
let CVID=0; const canvasJobs=[];   // large-N canvas draw jobs, rebuilt each render
function scatter(cfg){
  const W=cfg.W||520,H=cfg.H||360,ml=64,mr=22,mt=16,mb=48,pw=W-ml-mr,ph=H-mt-mb;
  const rows=cfg.rows.filter(g=>{const x=cfg.x(g),y=cfg.y(g);return x!=null&&isFinite(x)&&y!=null&&isFinite(y);});
  if(!rows.length) return '<p class="muted">no data for this panel</p>';
  let xs=rows.map(cfg.x), ys=rows.map(cfg.y);
  // Loop, not Math.min(...xs): spreading a huge array as arguments overflows the
  // call stack on whole-genome runs (tens of thousands of points).
  let xmin=Infinity,xmax=-Infinity,ymin=Infinity,ymax=-Infinity;
  for(let k=0;k<xs.length;k++){ const xv=xs[k],yv=ys[k];
    if(xv<xmin)xmin=xv; if(xv>xmax)xmax=xv; if(yv<ymin)ymin=yv; if(yv>ymax)ymax=yv; }
  (cfg.refs||[]).forEach(r=>{ if(r.x!=null&&isFinite(r.x)){xmin=Math.min(xmin,r.x);xmax=Math.max(xmax,r.x);} if(r.y!=null&&isFinite(r.y)){ymin=Math.min(ymin,r.y);ymax=Math.max(ymax,r.y);} });
  if(cfg.y0) ymin=Math.min(ymin,0);
  const xr=(xmax-xmin)||1, yr=(ymax-ymin)||1; xmin-=xr*0.04;xmax+=xr*0.04;ymin-=yr*0.04;ymax+=yr*0.08;
  const X=v=>ml+(v-xmin)/(xmax-xmin)*pw, Y=v=>mt+ph*(1-(v-ymin)/(ymax-ymin));
  let s=`<svg viewBox="0 0 ${W} ${H}">`;
  for(let i=0;i<=4;i++){ const yv=ymin+(ymax-ymin)*i/4, xv=xmin+(xmax-xmin)*i/4;
    s+=`<line x1="${ml}" y1="${Y(yv).toFixed(1)}" x2="${ml+pw}" y2="${Y(yv).toFixed(1)}" stroke="var(--grid)" stroke-width="0.5"/>`;
    s+=`<text x="${ml-8}" y="${Y(yv).toFixed(1)}" font-size="10" fill="var(--muted)" text-anchor="end" dominant-baseline="middle">${(cfg.yfmt||(v=>v.toFixed(2)))(yv)}</text>`;
    s+=`<text x="${X(xv).toFixed(1)}" y="${mt+ph+16}" font-size="10" fill="var(--muted)" text-anchor="middle">${(cfg.xfmt||(v=>v.toFixed(2)))(xv)}</text>`; }
  (cfg.refs||[]).forEach(r=>{ const c=r.c||"var(--line)";
    if(r.x!=null&&!isFinite(r.x)) return; if(r.y!=null&&!isFinite(r.y)) return;
    if(r.diag){ s+=`<line x1="${X(Math.max(xmin,ymin)).toFixed(1)}" y1="${Y(Math.max(xmin,ymin)).toFixed(1)}" x2="${X(Math.min(xmax,ymax)).toFixed(1)}" y2="${Y(Math.min(xmax,ymax)).toFixed(1)}" stroke="${c}" stroke-width="1" stroke-dasharray="6,3"/>`; }
    else if(r.x!=null){ s+=`<line x1="${X(r.x).toFixed(1)}" y1="${mt}" x2="${X(r.x).toFixed(1)}" y2="${mt+ph}" stroke="${c}" stroke-width="1" stroke-dasharray="6,3"/>`; if(r.label) s+=`<text x="${X(r.x).toFixed(1)}" y="${mt+10}" font-size="9" fill="${c}" text-anchor="middle">${r.label}</text>`; }
    else if(r.y!=null){ s+=`<line x1="${ml}" y1="${Y(r.y).toFixed(1)}" x2="${ml+pw}" y2="${Y(r.y).toFixed(1)}" stroke="${c}" stroke-width="1" stroke-dasharray="6,3"/>`; if(r.label) s+=`<text x="${ml+pw}" y="${(Y(r.y)-4).toFixed(1)}" font-size="9" fill="${c}" text-anchor="end">${r.label}</text>`; } });
  // For large gene sets, draw the points on a canvas overlay instead of thousands
  // of SVG nodes (keeps the page responsive on whole-genome runs).
  const useCanvas = rows.length > 1200;
  // Optional 95% CI whiskers on the y-axis (SVG mode only; skipped for huge N).
  if(cfg.yerr && !useCanvas){ rows.forEach(g=>{ const e=cfg.yerr(g); if(!e) return;
    const lo=Math.max(ymin,Math.min(ymax,e[0])), hi=Math.max(ymin,Math.min(ymax,e[1])); const x=X(cfg.x(g));
    s+=`<line x1="${x.toFixed(1)}" y1="${Y(lo).toFixed(1)}" x2="${x.toFixed(1)}" y2="${Y(hi).toFixed(1)}" stroke="${cfg.color(g)}" stroke-width="1" opacity="0.35"/>`; }); }
  if(!useCanvas){ rows.forEach(g=>{ const r=cfg.size?cfg.size(g):3.2, sel=(g._i===selected);
    const shape=(cvd&&cfg.shape)?cfg.shape(g):"dot";
    s+=marker(X(cfg.x(g)),Y(cfg.y(g)),r,shape,cfg.color(g),sel,g._i,cfg.tip(g)); }); }
  s+=`<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt+ph}" stroke="var(--axis)" stroke-width="1.2"/><line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="var(--axis)" stroke-width="1.2"/>`;
  s+=`<text x="${ml+pw/2}" y="${H-6}" font-size="11" fill="var(--muted)" text-anchor="middle">${cfg.xlabel}</text>`;
  s+=`<text x="14" y="${mt+ph/2}" font-size="11" fill="var(--muted)" text-anchor="middle" transform="rotate(-90,14,${mt+ph/2})">${cfg.ylabel}</text></svg>`;
  if(!useCanvas) return s;
  const id=++CVID;
  const pts=rows.map(g=>({x:X(cfg.x(g)),y:Y(cfg.y(g)),c:cfg.color(g),i:g._i,tip:cfg.tip(g),r:cfg.size?cfg.size(g):3.2,sel:g._i===selected,sh:(cvd&&cfg.shape)?cfg.shape(g):"dot"}));
  canvasJobs.push({id,W,H,pts});
  const svgAbs = s.replace('<svg viewBox','<svg style="position:absolute;inset:0;width:100%;height:100%" viewBox');
  return `<div class="cvwrap" style="position:relative">${svgAbs}<canvas id="cv${id}" width="${W}" height="${H}" style="position:relative;width:100%;height:auto;display:block" aria-label="scatter of ${rows.length} genes; use the table for details"></canvas></div>`;
}
// Draw all registered canvas jobs (large-N scatters) with nearest-point hover/click.
function drawCanvases(){
  const root=getComputedStyle(document.documentElement);
  const rc=c=>{ const m=/var\((--[\w-]+)\)/.exec(c||""); return m?(root.getPropertyValue(m[1]).trim()||c):c; };
  canvasJobs.forEach(job=>{ const cv=document.getElementById("cv"+job.id); if(!cv) return;
    const ctx=cv.getContext("2d"); ctx.clearRect(0,0,job.W,job.H);
    const sel=rc("var(--sel)");
    job.pts.forEach(p=>{ ctx.globalAlpha=p.sel?1:0.6; ctx.fillStyle=rc(p.c); const rr=p.sel?p.r+2:p.r;
      ctx.beginPath();
      if(p.sh==="up"){ const h=rr*1.5; ctx.moveTo(p.x,p.y-h); ctx.lineTo(p.x-h,p.y+h*0.66); ctx.lineTo(p.x+h,p.y+h*0.66); ctx.closePath(); }
      else if(p.sh==="down"){ const h=rr*1.5; ctx.moveTo(p.x,p.y+h); ctx.lineTo(p.x-h,p.y-h*0.66); ctx.lineTo(p.x+h,p.y-h*0.66); ctx.closePath(); }
      else { ctx.arc(p.x,p.y,rr,0,6.2832); }
      ctx.fill(); if(p.sel){ ctx.globalAlpha=1; ctx.strokeStyle=sel; ctx.lineWidth=2; ctx.stroke(); } });
    ctx.globalAlpha=1;
    const nearest=(e,maxd)=>{ const r=cv.getBoundingClientRect(); const mx=(e.clientX-r.left)*job.W/r.width, my=(e.clientY-r.top)*job.H/r.height;
      let best=null,bd=maxd; for(const p of job.pts){ const d=(p.x-mx)*(p.x-mx)+(p.y-my)*(p.y-my); if(d<bd){bd=d;best=p;} } return best; };
    // stopPropagation so the bubbling event does not reach the #panels mousemove
    // delegate, whose else-branch would immediately hide the tooltip we just set.
    cv.onmousemove=e=>{ e.stopPropagation(); const p=nearest(e,120); if(p){ tip.innerHTML=p.tip; tip.style.opacity=1; tip.style.left=(e.clientX+12)+"px"; tip.style.top=(e.clientY+12)+"px"; } else tip.style.opacity=0; };
    cv.onmouseleave=()=>tip.style.opacity=0;
    // Only swallow the click when a point is actually hit; otherwise let it bubble
    // so an open help popover still closes via the document outside-click handler.
    cv.onclick=e=>{ const p=nearest(e,80); if(p&&genes.some(g=>g._i===p.i)){ e.stopPropagation(); selectGene(p.i); } };
  });
}
// A single mark: circle by default, or a direction triangle in CVD mode so that
// direction is encoded by shape as well as colour. Keeps class="mark" + data-*
// so the existing hover/click delegation and selection ring work unchanged.
function marker(cx,cy,r,shape,fill,sel,gi,tip){
  const rr=sel?r+2:r;
  const a=`class="mark" data-i="${gi}" data-tip="${tip}" fill="${fill}" opacity="${sel?1:0.65}"${sel?' stroke="var(--sel)" stroke-width="2"':''}`;
  if(shape==="up"){ const h=rr*1.5; return `<path ${a} d="M${cx.toFixed(1)},${(cy-h).toFixed(1)} L${(cx-h).toFixed(1)},${(cy+h*0.66).toFixed(1)} L${(cx+h).toFixed(1)},${(cy+h*0.66).toFixed(1)} Z"/>`; }
  if(shape==="down"){ const h=rr*1.5; return `<path ${a} d="M${cx.toFixed(1)},${(cy+h).toFixed(1)} L${(cx-h).toFixed(1)},${(cy-h*0.66).toFixed(1)} L${(cx+h).toFixed(1)},${(cy-h*0.66).toFixed(1)} Z"/>`; }
  return `<circle ${a} cx="${cx.toFixed(1)}" cy="${cy.toFixed(1)}" r="${rr.toFixed(1)}"/>`;
}
const rSize = g => Math.min(9,Math.max(2.6,Math.sqrt((g.total||1))*1.1));
// Unified colour rule across every panel: red = significant diversifying (pN/pS>1),
// blue = significant purifying (<1), grey = not significant. (--sig and --pos share a
// hex, so we must never colour "significance" red independently of direction.)
const dirColor = g => isSig(g) ? (isUp(g) ? "var(--pos)" : "var(--accent)") : "var(--ns)";
const sigColor = dirColor;
const ciTip = g => (g.ratioLo!=null&&isFinite(g.ratioLo)) ? ` [${fmt(g.ratioLo,2)}–${isFinite(g.ratioHi)?fmt(g.ratioHi,2):'∞'}]` : '';
// HTML-escape data-derived strings (gene/chrom names) before putting them in
// innerHTML / SVG text / attributes, so a name with < > & " can't break markup.
const hesc = s => String(s==null?"":s).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;");
const baseTip = g => `<b>${hesc(g.name)}</b> (${hesc(g.chrom)}:${g.start} ${g.strand})<br>pN/pS ${fmt(g.ratio,3)}${ciTip(g)} · p ${fmtP(g.p)} · ${stringency==='q'?(M.genomicControl?'q(GC)':'q'):'p(Bonf)'} ${fmtP(sigVal(g))}<br>${g.nonsyn}N / ${g.syn}S of ${g.total} SNPs${quar(g)?' · ⚠ repetitive':''}`;
const log2c = (v,lo,hi) => v==null||!isFinite(v)?null : Math.max(lo,Math.min(hi, Math.log2(v<=0?1e-3:v)));

// ── Manhattan ─────────────────────────────────────────────────────────
function panelManhattan(){
  const thr = pThreshold();
  const yv = g => metric==="ratio" ? g.ratio : (metric==="z" ? g.z : -Math.log10(Math.max(pStat(g),1e-300)));
  const rows = genes.filter(g=> metric==="ratio" ? (g.ratio!=null&&isFinite(g.ratio)&&g.total>0)
    : metric==="z" ? (g.z!=null&&isFinite(g.z)) : (pStat(g)!=null&&isFinite(pStat(g))));
  const refs = metric==="ratio" ? [{y:1,label:"pN/pS = 1",c:"var(--muted)"}]
    : metric==="z" ? [{y:0,label:"z = 0",c:"var(--muted)"},{y:1.96,label:"+1.96",c:"var(--line)"},{y:-1.96,label:"−1.96",c:"var(--line)"}]
    : (thr!=null?[{y:-Math.log10(Math.max(thr,1e-300)),label:(stringency==="q"?"BH":"Bonf")+" "+M.fdr,c:"var(--line)"}]:[]);
  const tog=(v,t)=>`<button class="tog ${metric===v?'on':''}" data-act="metric" data-v="${v}">${t}</button>`;
  return {title:"Manhattan", help:"manhattan", span:true,
    toolbar: tog("neglogp","−log10(p)")+tog("ratio","pN/pS")+tog("z","z(N)"),
    legend:`<span><i style="background:var(--pos)"></i>sig. diversifying (pN/pS&gt;1)</span><span><i style="background:var(--accent)"></i>sig. purifying (&lt;1)</span><span><i style="background:var(--ns)"></i>not significant</span>${metric==="z"?'<span>· z = standardized nonsyn excess</span>':''}`,
    svg: scatter({rows,
      W:900,x:g=>g.start,y:yv,xlabel:"genome position",ylabel:metric==="ratio"?"pN/pS":(metric==="z"?"z (nonsyn excess)":"−log10(p)"),
      xfmt:v=>Math.round(v), color:dirColor, shape:dirShape,
      size:rSize,tip:baseTip, y0:metric==="z", refs})};
}
// ── Volcano ───────────────────────────────────────────────────────────
function panelVolcano(){
  const thr=pThreshold();
  return {title:"Volcano — effect vs significance", help:"volcano",
    legend:`<span><i style="background:var(--pos)"></i>sig. diversifying (right)</span><span><i style="background:var(--accent)"></i>sig. purifying (left)</span><span><i style="background:var(--ns)"></i>not significant</span>`,
    svg: scatter({rows:genes.filter(g=>pStat(g)!=null&&isFinite(pStat(g))&&g.ratio!=null),
      x:g=>log2c(g.ratio,-6,6), y:g=>-Math.log10(Math.max(pStat(g),1e-300)),
      xlabel:"log2(pN/pS)  ←purifying · positive→", ylabel:M.genomicControl?"−log10(p, GC)":"−log10(p)",
      color:sigColor,shape:dirShape,size:rSize,tip:baseTip,y0:true,
      refs:[{x:0,label:"pN/pS=1",c:"var(--muted)"}].concat(thr!=null?[{y:-Math.log10(Math.max(thr,1e-300)),label:(stringency==="q"?"BH":"Bonf"),c:"var(--line)"}]:[])})};
}
// ── Power funnel ──────────────────────────────────────────────────────
function panelFunnel(){
  return {title:"Power funnel — pN/pS vs SNP count", help:"funnel",
    legend:`<span><i style="background:var(--pos)"></i>sig. diversifying</span><span><i style="background:var(--accent)"></i>sig. purifying</span><span><i style="background:var(--ns)"></i>not significant</span><span>· low-count genes scatter widely</span>`,
    svg: scatter({rows:genes.filter(g=>g.total>0&&g.ratio!=null&&isFinite(g.ratio)),
      x:g=>Math.log10(g.total),y:g=>g.ratio,xlabel:"total SNPs (log10)",ylabel:"pN/pS",
      xfmt:v=>Math.round(Math.pow(10,v)),color:dirColor,shape:dirShape,size:rSize,tip:baseTip,y0:true,
      yerr:g=>(g.ratioLo!=null&&isFinite(g.ratioLo))?[g.ratioLo, isFinite(g.ratioHi)?g.ratioHi:g.ratio*4]:null,
      refs:[{y:1,label:"pN/pS=1",c:"var(--muted)"},{y:S.gwRatio,label:"pooled",c:"var(--line)"}].concat(M.minSnps>1?[{x:Math.log10(M.minSnps),label:"min-snps",c:"var(--muted)"}]:[])})};
}
// ── Observed vs expected nonsyn fraction ──────────────────────────────
function panelObsExp(){
  return {title:"Observed vs expected N-fraction", help:"obsexp",
    legend:`<span><i style="background:var(--pos)"></i>sig. above (diversifying)</span><span><i style="background:var(--accent)"></i>sig. below (purifying)</span><span><i style="background:var(--ns)"></i>not significant</span>`,
    svg: scatter({rows:genes.filter(g=>g.total>0&&g.expN!=null&&isFinite(g.expN)),
      x:g=>g.expN,y:g=>g.nonsyn/g.total,xlabel:"expected N-fraction  N/(N+S)",ylabel:"observed nonsyn / total",
      color:g=>{ if(!isSig(g))return"var(--ns)"; return (g.nonsyn/g.total)>g.expN?"var(--pos)":"var(--accent)"; },
      shape:g=>!isSig(g)?"dot":((g.nonsyn/g.total)>g.expN?"up":"down"),
      size:rSize,tip:g=>baseTip(g)+`<br>exp N-frac ${fmt(g.expN,3)} · obs ${fmt(g.nonsyn/g.total,3)}`,
      refs:[{diag:true,label:"neutral",c:"var(--muted)"}]})};
}
// ── MK volcano (if --mk) ──────────────────────────────────────────────
function panelMK(){
  return {title:"McDonald-Kreitman — α vs significance", help:"mk",
    legend:`<span><i style="background:var(--pos)"></i>adaptive (α&gt;0)</span><span><i style="background:var(--accent)"></i>constrained (α&lt;0)</span>`,
    svg: scatter({rows:genes.filter(g=>g.alpha!=null&&isFinite(g.alpha)&&g.fisherP!=null&&isFinite(g.fisherP)),
      x:g=>g.alpha,y:g=>-Math.log10(Math.max(g.fisherP,1e-300)),xlabel:"α (proportion adaptive)",ylabel:"−log10(Fisher p)",
      color:g=>g.alpha>0?"var(--pos)":"var(--accent)",shape:g=>g.alpha>0?"up":"down",size:g=>Math.min(9,Math.max(2.6,Math.sqrt((g.dn+g.ds+g.pnMk+g.psMk)||1)*1.3)),
      tip:g=>`<b>${hesc(g.name)}</b><br>Dn ${g.dn} Ds ${g.ds} · Pn ${g.pnMk} Ps ${g.psMk}<br>NI ${fmt(g.ni,3)} · α ${fmt(g.alpha,3)} · Fisher p ${fmtP(g.fisherP)}`,y0:true,
      refs:[{x:0,label:"α=0",c:"var(--line)"}]})};
}
// ── pN/pS distribution ────────────────────────────────────────────────
function panelDist(){
  const vals=genes.map(g=>g.ratio).filter(v=>v!=null&&isFinite(v));
  if(!vals.length) return null;
  const edges=[0,0.2,0.4,0.6,0.8,1.0,1.5,2.0,1e9], labels=["<.2",".2-.4",".4-.6",".6-.8",".8-1","1-1.5","1.5-2","≥2"];
  const counts=labels.map((_,k)=>vals.filter(v=>v>=edges[k]&&v<edges[k+1]).length);
  const cmax=Math.max(1,...counts),W=520,H=300,ml=54,mr=18,mt=16,mb=54,pw=W-ml-mr,ph=H-mt-mb,bw=pw/counts.length;
  let s=`<svg viewBox="0 0 ${W} ${H}">`;
  for(let i=0;i<=4;i++){ const v=cmax*i/4,y=mt+ph*(1-i/4); s+=`<line x1="${ml}" y1="${y}" x2="${ml+pw}" y2="${y}" stroke="var(--grid)" stroke-width="0.5"/><text x="${ml-8}" y="${y}" font-size="10" fill="var(--muted)" text-anchor="end" dominant-baseline="middle">${Math.round(v)}</text>`; }
  counts.forEach((c,i)=>{ const x=ml+i*bw+3,bh=ph*c/cmax,y=mt+ph-bh,
    col=edges[i]>=1?(edges[i]<1.5?"var(--ns)":"var(--pos)"):"var(--accent)";   // bin straddling pN/pS=1 shown neutral
    s+=`<rect x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${(bw-6).toFixed(1)}" height="${bh.toFixed(1)}" rx="3" fill="${col}" opacity="0.85"/>`;
    s+=`<text x="${(x+(bw-6)/2).toFixed(1)}" y="${mt+ph+16}" font-size="9" fill="var(--muted)" text-anchor="middle">${labels[i]}</text>`; });
  s+=`<text x="${ml+pw/2}" y="${H-6}" font-size="11" fill="var(--muted)" text-anchor="middle">pN/pS</text></svg>`;
  return {title:"pN/pS distribution", help:"dist", legend:"", svg:s};
}
// ── Site frequency spectrum: pN/pS by allele frequency ────────────────
HELP.sfs={title:"pN/pS by allele frequency (SFS)",lead:"Splits pN/pS by how common each variant is; purifying selection keeps deleterious nonsynonymous variants rare, so the ratio falls as frequency rises.",
  x:"alt allele frequency bin (rare → common)",y:"pooled pN/pS of the SNPs in that frequency bin (1 = neutral)",
  read:["A downward slope left→right: purifying selection removing nonsynonymous variants before they rise in frequency.","A flat profile near the genome-wide value: little frequency-dependent selection.","A high first bar with low later bars is the classic deleterious-variants-stay-rare signature."],
  watch:["Needs a multi-sample cohort with a real frequency spread; with one sample every variant sits at AF 1 and the panel is uninformative.","Bins with very few SNPs (small n) are noisy — read the n under each bar."]};
function panelSFS(){
  const nn=S.sfsNonsyn||[], sy=S.sfsSyn||[], edges=M.sfsEdges||[];
  if(!nn.length||!edges.length) return null;
  const nonempty=nn.map((v,i)=>v+(sy[i]||0)).filter(t=>t>0).length;
  const label=(lo,hi)=> lo===0?"≤"+fmt(hi,2):fmt(lo,2)+"–"+fmt(hi,2);
  if(nonempty<2){ return {title:"pN/pS by allele frequency (SFS)", help:"sfs", legend:"",
    extra:`<p class="muted">All variants fall in a single allele-frequency bin (single sample or invariant AF), so the frequency spectrum is uninformative here. This panel becomes meaningful with a multi-sample cohort.</p>`, svg:""}; }
  const scale=(S.gwSSites>0&&S.gwNSites>0)?(S.gwSSites/S.gwNSites):1;
  const bins=nn.map((_,i)=>{ const lo=i?edges[i-1]:0, hi=edges[i], n=nn[i], s=sy[i]||0;
    // n>0,s==0 → all-nonsynonymous bin → +∞; n==0 → pN/pS 0 (no nonsyn); else finite.
    return {lo,hi,n,s,ratio:(s>0&&n>0)?(n/s)*scale:(n>0?Infinity:0)}; });
  const finite=bins.map(b=>b.ratio).filter(v=>isFinite(v)); const vmax=Math.max(1.2,...finite);
  const W=560,H=330,ml=54,mr=18,mt=16,mb=70,pw=W-ml-mr,ph=H-mt-mb,bw=pw/bins.length;
  const Y=v=>mt+ph*(1-Math.min(v,vmax)/vmax);
  let s=`<svg viewBox="0 0 ${W} ${H}">`;
  for(let i=0;i<=4;i++){ const v=vmax*i/4,y=mt+ph*(1-i/4); s+=`<line x1="${ml}" y1="${y.toFixed(1)}" x2="${ml+pw}" y2="${y.toFixed(1)}" stroke="var(--grid)" stroke-width="0.5"/><text x="${ml-8}" y="${y.toFixed(1)}" font-size="10" fill="var(--muted)" text-anchor="end" dominant-baseline="middle">${fmt(v,2)}</text>`; }
  s+=`<line x1="${ml}" y1="${Y(1).toFixed(1)}" x2="${ml+pw}" y2="${Y(1).toFixed(1)}" stroke="var(--muted)" stroke-width="1" stroke-dasharray="6,3"/><text x="${(ml+pw).toFixed(1)}" y="${(Y(1)-4).toFixed(1)}" font-size="9" fill="var(--muted)" text-anchor="end">neutral 1</text>`;
  bins.forEach((b,i)=>{ const x=ml+i*bw+4, bh=ph*Math.min(isFinite(b.ratio)?b.ratio:vmax,vmax)/vmax, y=mt+ph-bh, col=b.ratio>1?"var(--pos)":"var(--accent)";
    s+=`<rect class="mark" data-tip="AF ${label(b.lo,b.hi)}<br>pN/pS ${isFinite(b.ratio)?fmt(b.ratio,3):"∞"}<br>${b.n}N / ${b.s}S SNPs" x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${(bw-8).toFixed(1)}" height="${bh.toFixed(1)}" rx="3" fill="${col}" opacity="0.85"/>`;
    s+=`<text x="${(x+(bw-8)/2).toFixed(1)}" y="${mt+ph+15}" font-size="9" fill="var(--muted)" text-anchor="middle">${label(b.lo,b.hi)}</text>`;
    s+=`<text x="${(x+(bw-8)/2).toFixed(1)}" y="${mt+ph+29}" font-size="8" fill="var(--muted)" text-anchor="middle">n=${b.n+b.s}</text>`; });
  s+=`<text x="${(ml+pw/2).toFixed(1)}" y="${H-8}" font-size="11" fill="var(--muted)" text-anchor="middle">alt allele frequency →</text></svg>`;
  return {title:"pN/pS by allele frequency (SFS)", help:"sfs",
    legend:`<span>purifying selection makes pN/pS fall as allele frequency rises</span>`, svg:s};
}

// ── Polymorphism vs divergence reconciliation (--divergence) ──────────
const hasDiv = genes.some(g=>g.div!=null&&isFinite(g.div));
function panelRecon(){
  const rows=genes.filter(g=>g.div!=null&&isFinite(g.div)&&g.ratio!=null&&isFinite(g.ratio));
  const pastPos=rows.filter(g=>g.div>1&&g.ratio<1).sort((a,b)=>b.div-a.div);
  const cand = pastPos.length
    ? `<div class="candlist"><b>Past-positive candidates</b> (dN/dS&gt;1, pN/pS&lt;1) — `+
      pastPos.slice(0,15).map(g=>`<span class="cand" data-i="${g._i}">${hesc(g.name)} <small>${fmt(g.div,2)}/${fmt(g.ratio,2)}</small></span>`).join("")+
      (pastPos.length>15?`<span class="muted">+${pastPos.length-15} more</span>`:"")+`</div>`
    : `<div class="candlist muted">No past-positive candidates in this set.</div>`;
  return {title:`Polymorphism vs divergence (${rows.length} genes matched)`, help:"recon", span:true,
    legend:`<span><i style="background:var(--pos)"></i>past positive (fixed)</span><span><i style="background:var(--s3)"></i>diversifying</span><span><i style="background:var(--s7)"></i>relaxed/recent</span><span><i style="background:var(--accent)"></i>purifying</span>`,
    extra: cand,
    svg: scatter({rows, x:g=>g.ratio, y:g=>g.div, xlabel:"pN/pS (within-sample polymorphism)", ylabel:"dN/dS (divergence)",
      color:g=>{ const x=g.ratio,y=g.div; if(y>1&&x<1)return"var(--pos)"; if(x>1&&y>1)return"var(--s3)"; if(x>1&&y<1)return"var(--s7)"; return"var(--accent)"; },
      size:g=>{ const base=rSize(g); return (g.div>1&&g.ratio<1)?base+1.5:base; },
      tip:g=>`<b>${hesc(g.name)}</b><br>pN/pS ${fmt(g.ratio,3)} · dN/dS ${fmt(g.div,3)}<br>${g.nonsyn}N/${g.syn}S SNPs`,
      refs:[{diag:true,label:"concordance",c:"var(--line)"},{x:1,label:"pN/pS=1",c:"var(--muted)"},{y:1,label:"dN/dS=1",c:"var(--muted)"}]})};
}
// ── Selection-regime census (click a regime to filter the table) ──────
function panelCensus(){
  const counts={positive:0,purifying:0,ns:0}; genes.forEach(g=>counts[regime(g)]++);
  const chip=(k,n)=>`<button class="chip regime ${regimeFilter===k?'on':''}" data-regime="${k}"><i style="background:${REGIMES[k]}"></i>${RLAB[k]} <b>${n}</b></button>`;
  const all=`<button class="chip regime ${regimeFilter==null?'on':''}" data-regime="">all <b>${genes.length}</b></button>`;
  return {title:"Selection regimes", help:"census", span:true,
    legend:`<span>click a regime to filter the table below · "not significant" = insufficient SNPs to reject neutrality (may still be selected)</span>`,
    extra:`<div class="census">${all}${chip("positive",counts.positive)}${chip("purifying",counts.purifying)}${chip("ns",counts.ns)}</div>`,
    svg:""};
}

// ── Significant-hits shortlist (always shown; explicit empty state) ────
function panelHits(){
  const hits=genes.filter(isSig).slice().sort((a,b)=>{ const sa=sigVal(a),sb=sigVal(b);
    const na=(sa==null||!isFinite(sa)),nb=(sb==null||!isFinite(sb)); if(na&&nb)return 0; if(na)return 1; if(nb)return -1; return sa-sb; });
  const body = hits.length
    ? `<div class="candlist">`+hits.map(g=>`<span class="cand" data-i="${g._i}">${isUp(g)?"▲":"▼"} ${hesc(g.name)} <small>${fmt(g.ratio,2)} · ${stringency==="q"?"q":"pB"} ${fmtP(sigVal(g))}${quar(g)?" ⚠":""}</small></span>`).join("")+`</div>`
    : `<div class="candlist muted">No gene rejects neutrality at ${stringency==="q"?"BH-FDR":"Bonferroni"} &lt; ${M.fdr}. In clonal M. tuberculosis this is common — most genes are underpowered rather than neutral.</div>`;
  return {title:`Significant hits (${hits.length})`, help:"hits", span:true, legend:"", extra:body, svg:""};
}
// ── P-value QQ + genomic inflation λ (λ comes from the backend) ───────
function panelQQ(){
  // Use the log-space −log10(p) so genes whose exact p underflowed to 0 still plot.
  const rows=genes.filter(g=>g.p!=null&&isFinite(g.p))
    .map(g=>({p:g.p, nlp:(g.nlp!=null&&isFinite(g.nlp))?g.nlp:-Math.log10(Math.max(g.p,1e-300))}))
    .sort((a,b)=>a.p-b.p);
  const m=rows.length; if(m<2) return null;
  rows.forEach((r,i)=>{ r.exp=-Math.log10((i+0.5)/m); r.obs=r.nlp; r._i=-1; });
  const lambda=(M.lambda!=null&&isFinite(M.lambda))?M.lambda:null;
  return {title:"P-value QQ"+(lambda!=null?` · λ=${fmt(lambda,2)}`:""), help:"qq",
    legend:`<span>points above the y=x line = more significant genes than a uniform null${lambda!=null?` · λ (inflation) = ${fmt(lambda,2)}`:""}</span>`,
    svg: scatter({rows, x:r=>r.exp, y:r=>r.obs, xlabel:"expected −log10(p)  (uniform null)", ylabel:"observed −log10(p)",
      color:()=>"var(--accent)", size:()=>3.2, tip:r=>`p ${fmtP(r.p)}<br>exp ${fmt(r.exp,2)} · obs ${fmt(r.obs,2)}`,
      refs:[{diag:true,label:"y = x",c:"var(--muted)"}]})};
}
// ── Ranked lollipop of top genes by selection signal ──────────────────
function panelLollipop(){
  const tested=genes.filter(g=>g.ratio!=null&&isFinite(g.ratio)&&g.total>0);
  if(!tested.length) return null;
  const absl=g=>Math.abs(Math.log2(g.ratio<=0?1e-3:g.ratio));
  const ranked=tested.slice().sort((a,b)=>{ const sa=sigVal(a),sb=sigVal(b);
    const na=(sa==null||!isFinite(sa)),nb=(sb==null||!isFinite(sb));
    if(!na&&!nb&&sa!==sb) return sa-sb; if(na!==nb) return na?1:-1; return absl(b)-absl(a); }).slice(0,18);
  const n=ranked.length,W=900,rowH=22,mt=8,mb=42,ml=160,mr=104,H=mt+mb+n*rowH,pw=W-ml-mr,lo=-4,hi=4;
  const X=v=>ml+(Math.max(lo,Math.min(hi,v))-lo)/(hi-lo)*pw, x0=X(0);
  let s=`<svg viewBox="0 0 ${W} ${H}">`;
  [-2,-1,0,1,2].forEach(t=>{ const x=X(t); s+=`<line x1="${x.toFixed(1)}" y1="${mt}" x2="${x.toFixed(1)}" y2="${mt+n*rowH}" stroke="${t===0?'var(--muted)':'var(--grid)'}" stroke-width="${t===0?1:0.5}"${t===0?' stroke-dasharray="4,3"':''}/>`;
    s+=`<text x="${x.toFixed(1)}" y="${mt+n*rowH+15}" font-size="9" fill="var(--muted)" text-anchor="middle">${Math.pow(2,t).toFixed(t<0?2:0)}</text>`; });
  ranked.forEach((g,i)=>{ const y=mt+i*rowH+rowH/2, lv=Math.log2(g.ratio<=0?1e-3:g.ratio), x=X(lv), col=REGIMES[regime(g)], sel=(g._i===selected);
    s+=`<line x1="${x0.toFixed(1)}" y1="${y}" x2="${x.toFixed(1)}" y2="${y}" stroke="${col}" stroke-width="2" opacity="0.45"/>`;
    s+=`<circle class="mark" data-i="${g._i}" data-tip="${baseTip(g)}" cx="${x.toFixed(1)}" cy="${y}" r="${sel?6:4.5}" fill="${col}"${sel?' stroke="var(--sel)" stroke-width="2"':''}/>`;
    s+=`<text x="${ml-8}" y="${y}" font-size="10" fill="var(--fg)" text-anchor="end" dominant-baseline="middle">${hesc(g.name)}${quar(g)?" ⚠":""}</text>`;
    s+=`<text x="${(ml+pw+8).toFixed(1)}" y="${y}" font-size="9" fill="var(--muted)" dominant-baseline="middle">${fmt(g.ratio,2)} · ${isSig(g)?(stringency==="q"?"q":"pB")+" "+fmtP(sigVal(g)):"ns"}</text>`; });
  s+=`<text x="${(ml+pw/2).toFixed(1)}" y="${H-6}" font-size="11" fill="var(--muted)" text-anchor="middle">pN/pS (log2 scale) — stem anchored at neutral (1)</text></svg>`;
  return {title:"Top genes by selection signal", help:"lollipop", span:true,
    legend:`<span><i style="background:var(--pos)"></i>positive</span><span><i style="background:var(--accent)"></i>purifying</span><span><i style="background:var(--ns)"></i>not significant</span><span>⚠ repetitive/PE-PPE</span>`,
    svg:s};
}
HELP.hits={title:"Significant hits",lead:"A shortlist of exactly the genes that reject neutrality at the current stringency, ranked by evidence.",x:"",y:"",
  read:["Each chip is a significant gene: ▲ pN/pS&gt;1 (diversifying), ▼ &lt;1 (purifying), with its ratio and q / p(Bonferroni).","Click a chip to highlight that gene across every panel and the table.","Toggle FDR(BH)↔Bonferroni in the top bar to widen or tighten this list."],
  watch:["⚠ marks repetitive / PE-PPE genes whose SNPs are often mapping artefacts — treat their hits sceptically."]};
HELP.lollipop={title:"Top genes by selection signal",lead:"The strongest-signal genes ranked together, so the leading candidates are visible at a glance.",
  x:"pN/pS on a log2 scale (stem anchored at the neutral value 1); left = purifying, right = diversifying",y:"one gene per row, ordered by significance then by effect size",
  read:["Longest bars to the right = biggest nonsynonymous excess; to the left = strongest constraint.","Colour marks the regime; the q / p(Bonferroni) at the right tells you if the bar is significant.","Click any dot to select that gene everywhere."],
  watch:["Ranking favours well-powered genes; a long bar on very few SNPs can still be noise — check the count and q-value."]};

// ── Core vs repetitive pooled ratio (only when repeats exist) ─────────
HELP.corerep={title:"Core vs repetitive pN/pS",lead:"The genome-wide ratio split into core genes and repetitive ones (PE/PPE/PGRS, IS elements), which map poorly.",x:"",y:"",
  read:["A repetitive pN/pS far above core is the fingerprint of mapping/alignment artefacts, not real selection.","If they agree, the repeats are not distorting the genome-wide estimate."],
  watch:["--exclude-repetitive drops the repetitive set from the pooled estimate and the neutrality-test family."]};
function panelCoreRep(){
  const c=S.coreGw, r=S.repGw; if(!r) return null;
  const row=(lab,g)=> g?`<div><b>${lab}:</b> pN/pS <b>${fmt(g.ratio,3)}</b> <small style="color:var(--muted)">(${g.nonsyn}N / ${g.syn}S SNPs)</small></div>`:`<div><b>${lab}:</b> —</div>`;
  return {title:"Core vs repetitive pN/pS", help:"corerep", legend:"",
    extra:`<div style="display:flex;flex-direction:column;gap:5px;font-size:.9rem">${row("core",c)}${row("repetitive",r)}<div class="muted">${M.excludeRepetitive?"repetitive genes are excluded from the genome-wide estimate and the test":"rerun with --exclude-repetitive to drop repeats from the pooled estimate + test"}</div></div>`, svg:""};
}

// ── Render all panels ─────────────────────────────────────────────────
function renderPanels(){
  hideInfo();   // any open help popover is anchored to a button we are about to destroy
  if(typeof tip!=="undefined"&&tip) tip.style.opacity=0;   // clear any stale tooltip before rebuild
  canvasJobs.length=0;   // discard previous render's canvas jobs before rebuilding
  if(!S.genesWithSnps){ $("#panels").innerHTML=`<section class="panel wide"><p class="muted">No genes carried usable SNPs, so there is nothing to plot. Check the VCF/GFF contig names and filters.</p></section>`; buildToc(); return; }
  const P=[panelCensus(),panelHits(),panelManhattan(),panelVolcano(),panelQQ()];
  if(M.mk) P.push(panelMK());
  if(hasDiv) P.push(panelRecon());
  P.push(panelLollipop(),panelFunnel(),panelObsExp());
  const sfs=panelSFS(); if(sfs) P.push(sfs);
  const cr=panelCoreRep(); if(cr) P.push(cr);
  const d=panelDist(); if(d) P.push(d);
  $("#panels").innerHTML = P.filter(Boolean).map((p,i)=>`<section id="sec-${p.help||('p'+i)}" class="panel${p.span?' wide':''}"><h2><span class="htitle">${p.title}</span>${p.help&&HELP[p.help]?`<button class="info" data-help="${p.help}" title="How to read this" aria-label="How to read this panel">i</button>`:""}</h2>${p.toolbar?`<div class="toolbar">${p.toolbar}</div>`:""}${p.svg?`<div>${p.svg}</div>`:""}${p.legend?`<div class="legend">${p.legend}</div>`:""}${p.extra||""}</section>`).join("");
  drawCanvases();   // paint any large-N canvas panels now that they are in the DOM
  buildToc();       // the left table of contents tracks the (dynamic) panel list
}
// ── Left table of contents with scroll-spy ────────────────────────────
function buildToc(){
  const secs=[["verdict","Summary"]];
  document.querySelectorAll("#panels section.panel").forEach(sec=>{
    if(sec.id){ const t=sec.querySelector(".htitle"); secs.push([sec.id, t?t.textContent:sec.id]); } });
  secs.push(["sec-table","Per-gene table"]);
  $("#toc").innerHTML='<div class="toc-h">Contents</div>'+
    secs.map(([id,t])=>`<a href="#${id}" data-target="${id}" title="${t}">${t}</a>`).join("");
  updateTocActive();
}
// Mark the last section whose top has scrolled above a line near the viewport top.
function updateTocActive(){
  const links=[...document.querySelectorAll("#toc a")]; if(!links.length) return;
  const y=window.scrollY+90; let current=links[0];
  for(const a of links){ const el=document.getElementById(a.dataset.target);
    if(!el) continue; const top=el.getBoundingClientRect().top+window.scrollY;
    if(top<=y) current=a; else break; }
  links.forEach(a=>a.classList.toggle("active", a===current));
}
let tocTick=false;
window.addEventListener("scroll", ()=>{ if(tocTick) return; tocTick=true;
  requestAnimationFrame(()=>{ tocTick=false; updateTocActive(); }); }, {passive:true});
document.querySelector("#toc").addEventListener("click", e=>{ const a=e.target.closest("a"); if(!a) return;
  e.preventDefault(); const el=document.getElementById(a.dataset.target); if(el) el.scrollIntoView({behavior:"smooth",block:"start"}); });
function fixSpans(){}  // wide panels now use the .wide class

// ── Interaction: hover tooltip + click-to-select (delegated) ──────────
$("#panels").addEventListener("mousemove", e=>{ const t=e.target;
  if(t.classList&&t.classList.contains("mark")){ tip.innerHTML=t.dataset.tip; tip.style.opacity=1; tip.style.left=(e.clientX+12)+"px"; tip.style.top=(e.clientY+12)+"px"; } else tip.style.opacity=0; });
$("#panels").addEventListener("mouseleave", ()=>tip.style.opacity=0);
$("#panels").addEventListener("click", e=>{ const t=e.target;
  if(t.classList&&t.classList.contains("mark")){ const i=+t.dataset.i; if(genes.some(g=>g._i===i)) selectGene(i); return; }
  const cand=t.closest&&t.closest(".cand"); if(cand){ selectGene(+cand.dataset.i); return; }
  const chip=t.closest&&t.closest(".chip.regime"); if(chip){ regimeFilter=chip.dataset.regime||null; renderPanels(); fixSpans(); renderTable(); } });

// Bring the selected gene's row into view, even when the table is virtualized
// (the row may not be materialised, so scroll the container to it first).
function scrollToSelected(){
  if(selected==null) return;
  if(tableRows.length>VIRT_MIN){ const idx=tableRows.findIndex(g=>g._i===selected);
    if(idx>=0){ const wrap=$("#tblwrap"), y=idx*ROWH;
      if(y<wrap.scrollTop || y>wrap.scrollTop+wrap.clientHeight-ROWH){ wrap.scrollTop=Math.max(0,y-wrap.clientHeight/2); }
      renderTableWindow(); } }
  const row=document.querySelector(`#tbl tbody tr[data-i="${selected}"]`); if(row) row.scrollIntoView({block:"nearest"});
}
function selectGene(i){ selected = (selected===i)?null:i;
  // If a filter/search would hide the newly-selected gene, relax it so the
  // highlighted point always has a findable table row.
  if(selected!=null && !filtered().some(g=>g._i===selected)){ regimeFilter=null; $("#search").value=""; }
  renderPanels(); fixSpans(); renderTable();
  scrollToSelected(); }

// ── Table ─────────────────────────────────────────────────────────────
let cols=[["name","Gene"],["chrom","Chrom"],["start","Start"],["strand","±"],["length_bp","Len"],
  ["pn","pN"],["ps","pS"],["ratio",M.ratioName],["ci","95% CI"],["dir","dir"],["nonsyn","N"],["syn","S"],["total","SNPs"],
  ["expN","ExpN"],["z","z(N)"],["p","p"],["q","q(BH)"],["bonf","p(Bonf)"]];
if(M.genomicControl) cols=cols.concat([["pGc","p(GC)"],["qGc","q(GC)"]]);
if(M.mk) cols=cols.concat([["dn","Dn"],["ds","Ds"],["pnMk","Pn"],["psMk","Ps"],["ni","NI"],["alpha","α"],["dos","DoS"],["fisherP","Fisher"]]);
const pcols=new Set(["p","q","bonf","fisherP","pGc","qGc"]); const icols=new Set(["start","length_bp","nonsyn","syn","total","dn","ds","pnMk","psMk"]);
let sortKey="p",sortDesc=false;
$("#tbl thead").innerHTML="<tr>"+cols.map(c=>`<th data-k="${c[0]}" scope="col">${c[1]}</th>`).join("")+"</tr>";
function cellText(g,k){
  if(k==="dir"){ if(g.ratio!=null&&isFinite(g.ratio)) return g.ratio>1?"▲":(g.ratio<1?"▼":"·");
    return (g.total>0&&(g.nonsyn||0)!==(g.syn||0)) ? (g.nonsyn>g.syn?"▲":"▼") : "·"; }
  if(k==="ci"){ return (g.ratioLo!=null&&isFinite(g.ratioLo)) ? `${fmt(g.ratioLo,2)}–${isFinite(g.ratioHi)?fmt(g.ratioHi,2):"∞"}` : "NA"; }
  let v=g[k]; if(typeof v==="string") return k==="name"?hesc(v)+(quar(g)?'<span class="badge">rep</span>':""):hesc(v);
  if(icols.has(k)) return v==null?"NA":v;
  return pcols.has(k)?fmtP(v):fmt(v,(k==="ratio"||k==="ni"||k==="alpha"||k==="pn"||k==="ps"||k==="expN")?4:2);
}
let tableRows=[], ROWH=29; const VIRT_MIN=400;   // virtualize the table above VIRT_MIN rows
function rowHtml(g, idx){
  // Zebra by ABSOLUTE row index (idx) so striping is stable under virtualization,
  // where a leading spacer <tr> would otherwise flip :nth-child parity.
  return `<tr data-i="${g._i}" class="${isSig(g)?'sig ':''}${idx%2?'zebra ':''}${g._i===selected?'sel':''}">`+cols.map(c=>{
    const col = c[0]==="dir" ? (isUp(g)?"var(--pos)":"var(--accent)") : null;
    return `<td${col?` style="color:${col}"`:""}>${cellText(g,c[0])}</td>`; }).join("")+"</tr>";
}
function renderTableWindow(){
  const rows=tableRows, tb=$("#tbl tbody");
  if(rows.length<=VIRT_MIN){ tb.innerHTML=rows.map((g,i)=>rowHtml(g,i)).join(""); return; }
  // Only materialise the rows near the viewport; pad with spacer rows so the
  // scrollbar still reflects the full height.
  const wrap=$("#tblwrap"), st=wrap.scrollTop, vh=wrap.clientHeight||560;
  const first=Math.max(0,Math.floor(st/ROWH)-6), count=Math.ceil(vh/ROWH)+12;
  const slice=rows.slice(first,first+count);
  const padTop=first*ROWH, padBot=Math.max(0,(rows.length-first-slice.length)*ROWH);
  const nc=cols.length;
  tb.innerHTML=(padTop>0?`<tr class="pad" style="height:${padTop}px"><td colspan="${nc}"></td></tr>`:"")
    +slice.map((g,i)=>rowHtml(g,first+i)).join("")
    +(padBot>0?`<tr class="pad" style="height:${padBot}px"><td colspan="${nc}"></td></tr>`:"");
  // Calibrate ROWH from a real materialised row so spacer heights track the true
  // layout (matters for the scrollbar over thousands of rows).
  const dr=tb.querySelector("tr[data-i]"); if(dr){ const h=dr.getBoundingClientRect().height; if(h>4) ROWH=h; }
}
function renderTable(){
  let rows=filtered();
  const sk = (sortKey==="dir"||sortKey==="ci") ? "ratio" : sortKey;   // 'dir'/'ci' derive from ratio
  rows.sort((a,b)=>{ let x=a[sk],y=b[sk]; const xn=(x==null||(typeof x==="number"&&!isFinite(x))),yn=(y==null||(typeof y==="number"&&!isFinite(y)));
    if(xn&&yn)return 0; if(xn)return 1; if(yn)return -1;
    if(typeof x==="string")return sortDesc?y.localeCompare(x):x.localeCompare(y); return sortDesc?y-x:x-y; });
  tableRows=rows; renderTableWindow();
  document.querySelectorAll("#tbl thead th").forEach(th=>{ th.classList.toggle("sorted",th.dataset.k===sortKey); th.classList.toggle("desc",th.dataset.k===sortKey&&sortDesc); });
  $("#tableCount").textContent=`${rows.length} / ${genes.length} genes`+(rows.length>VIRT_MIN?" (virtualized)":"");
}
// Re-materialise the visible window as the table scrolls (virtualized mode only).
$("#tblwrap").addEventListener("scroll", ()=>{ if(tableRows.length>VIRT_MIN) renderTableWindow(); }, {passive:true});
$("#tbl thead").addEventListener("click",e=>{ const th=e.target.closest("th"); if(!th)return; const k=th.dataset.k; if(k===sortKey)sortDesc=!sortDesc; else{sortKey=k;sortDesc=false;} renderTable(); });
$("#tbl tbody").addEventListener("click",e=>{ const tr=e.target.closest("tr"); if(tr) selectGene(+tr.dataset.i); });
$("#search").addEventListener("input", renderTable);

// ── Panel toolbar (metric toggle) via delegation ──────────────────────
$("#panels").addEventListener("click", e=>{ const b=e.target.closest("button.tog"); if(!b)return;
  if(b.dataset.act==="metric"){ metric=b.dataset.v; renderPanels(); fixSpans(); } });

// ── Stringency toggle ─────────────────────────────────────────────────
$("#strTog").addEventListener("click", ()=>{ stringency = stringency==="q"?"bonf":"q";
  $("#strTog").textContent = stringency==="q"?"FDR (BH)":"Bonferroni"; renderVerdict(); renderCards(); renderPanels(); fixSpans(); renderTable(); });

// ── Theme toggle ──────────────────────────────────────────────────────
$("#themeTog").addEventListener("click", ()=>{ const r=document.documentElement;
  const cur=r.getAttribute("data-theme")||(matchMedia("(prefers-color-scheme:dark)").matches?"dark":"light");
  r.setAttribute("data-theme", cur==="dark"?"light":"dark");
  drawCanvases(); });   // canvas colours are baked at paint time — re-resolve them for the new theme

// ── Colour-blind (CVD) mode: Okabe-Ito palette + direction shapes ─────
$("#cvdTog").addEventListener("click", ()=>{ cvd=!cvd; const r=document.documentElement;
  if(cvd) r.setAttribute("data-cvd",""); else r.removeAttribute("data-cvd");
  const b=$("#cvdTog"); b.classList.toggle("on",cvd); b.setAttribute("aria-pressed",cvd?"true":"false");
  renderPanels(); });   // repaint so marks pick up their direction shapes

// ── Export CSV / JSON of the filtered view ────────────────────────────
function filtered(){ const f=($("#search").value||"").toLowerCase();
  return genes.filter(g=>(!f||g.name.toLowerCase().includes(f)||(g.chrom||"").toLowerCase().includes(f)) && (!regimeFilter||regime(g)===regimeFilter)); }
function dl(name,txt,type){ const b=new Blob([txt],{type}); const a=document.createElement("a"); a.href=URL.createObjectURL(b); a.download=name; a.click(); URL.revokeObjectURL(a.href); }
$("#expCsv").addEventListener("click", ()=>{ const keys=cols.map(c=>c[0]).filter(k=>k!=="dir"&&k!=="ci").concat(["ratioLo","ratioHi"]);
  const cell=v=>{ const s=(v==null||(typeof v==="number"&&!isFinite(v)))?"NA":String(v); return /[",\n]/.test(s)?'"'+s.replace(/"/g,'""')+'"':s; };
  const head=keys.map(cell).join(","); const body=filtered().map(g=>keys.map(k=>cell(g[k])).join(",")).join("\n");
  dl("eskaks_genes.csv", head+"\n"+body, "text/csv"); });
$("#expJson").addEventListener("click", ()=>{ dl("eskaks_genes.json", JSON.stringify(filtered().map(g=>{const o={...g};delete o._i;return o;}),null,1), "application/json"); });

// ── Keyboard navigation: ↑/↓ step the selected gene, Esc clears ───────
document.addEventListener("keydown", e=>{
  const ae=document.activeElement; if(ae&&(ae.tagName==="INPUT"||ae.tagName==="TEXTAREA")) return;
  if(e.key==="ArrowDown"||e.key==="ArrowUp"){
    const nsk=(sortKey==="dir"||sortKey==="ci")?"ratio":sortKey;   // match renderTable's derived-column mapping
    const list=filtered().slice().sort((a,b)=>{ let x=a[nsk],y=b[nsk];
      const xn=(x==null||(typeof x==="number"&&!isFinite(x))),yn=(y==null||(typeof y==="number"&&!isFinite(y)));
      if(xn&&yn)return 0; if(xn)return 1; if(yn)return -1; if(typeof x==="string")return sortDesc?y.localeCompare(x):x.localeCompare(y); return sortDesc?y-x:x-y; });
    if(!list.length) return; e.preventDefault();
    let idx=list.findIndex(g=>g._i===selected);
    idx = e.key==="ArrowDown" ? Math.min(list.length-1, idx+1) : Math.max(0, (idx<0?0:idx-1));
    selected=list[idx]._i; renderPanels(); renderTable();
    scrollToSelected();
  } else if(e.key==="Escape" && !openInfo && selected!=null){ selected=null; renderPanels(); renderTable(); }
});
// ── Print / Save-PDF (opens collapsibles, then the print dialog) ──────
$("#printBtn").addEventListener("click", ()=>{ hideInfo();
  document.querySelectorAll("details").forEach(d=>d.open=true); window.print(); });

renderVerdict(); renderGuide(); renderCards(); renderPanels(); fixSpans(); renderTable();
"##;

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
    html.push_str(FASTA_BODY);
    html.push_str("<script>\nconst DATA = ");
    html.push_str(&data);
    html.push_str(";\n");
    html.push_str(FASTA_SCRIPT);
    html.push_str("</script>\n</body>\n</html>\n");

    let mut file = File::create(&output_path)?;
    file.write_all(html.as_bytes())?;
    Ok(output_path)
}

const FASTA_BODY: &str = r#"<div class="wrap">
<h1>eskaks — dN/dS report</h1>
<div class="sub" id="meta"></div>
<div class="topbar">
  <button class="btn" id="expJson">⤓ JSON</button>
  <button class="btn" id="themeTog" title="Toggle light/dark">◑ Theme</button>
</div>
<div class="cards" id="cards"></div>
<div id="sections"></div>
<footer>Generated by eskaks · self-contained, no external assets</footer>
</div>
<div class="tip" id="tip"></div>
"#;

const FASTA_SCRIPT: &str = r##"
const $ = s => document.querySelector(s);
const M = DATA.meta, tip = $("#tip");
// HTML-escape genome/lineage/group names before putting them in SVG text, tooltips
// or legends, so a name with < > & " can't break markup or inject HTML.
const hesc = s => String(s==null?"":s).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;");
const fmt = (v,d=4) => v==null||!isFinite(v) ? "NA" : Number(v).toFixed(d);
// PathoGenOmics-Lab mycolorsTB — M. tuberculosis lineage colors (keyed by name).
const MYCOLORS = {A1:"#d1ae00",A2:"#8ef5c8",A3:"#73c2ff",A4:"#ff9cdb",L1:"#ff3091",L2:"#001aff",
  L3:"#8a0bd2",L4:"#ff0000",L5:"#995200",L6:"#1eb040",L7:"#fbff00",L8:"#ff9d00",L9:"#37ff30",L10:"#8fbda1"};
const SER = ["--s1","--s2","--s3","--s4","--s5","--s6","--s7","--s8"];

$("#meta").textContent = `model ${M.model} · ${M.totalPairs} pairs (${M.validPairs} valid)`;
const card = (k,v) => `<div class="card"><div class="k">${k}</div><div class="v">${v}</div></div>`;
$("#cards").innerHTML =
  card("Pairs", M.totalPairs) + card("Valid pairs", M.validPairs) +
  card("Pooled dN/dS", fmt(M.pooled)) + card("Mean dN / dS", fmt(M.meanDn,3)+" / "+fmt(M.meanDs,3));

const sec = $("#sections");
function addSection(title, id){ const d=document.createElement("section"); d.innerHTML=`<h2>${title}</h2><div style="overflow-x:auto" id="${id}"></div>`; sec.appendChild(d); return $("#"+id); }
function wireTip(c){
  c.addEventListener("mousemove", e=>{ const info=e.target.getAttribute && e.target.getAttribute("data-tip");
    if(info){ tip.innerHTML=info; tip.style.opacity=1; tip.style.left=(e.clientX+12)+"px"; tip.style.top=(e.clientY+12)+"px"; } else tip.style.opacity=0; });
  c.addEventListener("mouseleave", ()=>tip.style.opacity=0);
}
function grid(s,ml,pw,mt,ph,ymax,Y,dec){ for(let i=0;i<=5;i++){ const val=ymax*i/5,y=Y(val);
  s.v+=`<line x1="${ml}" y1="${y}" x2="${ml+pw}" y2="${y}" stroke="var(--border)" stroke-width="0.5"/><text x="${ml-8}" y="${y}" font-size="10" fill="var(--muted)" text-anchor="end" dominant-baseline="middle">${val.toFixed(dec)}</text>`; } }
function neutral(s,ml,pw,ymax,Y){ if(1<=ymax){ const y=Y(1);
  s.v+=`<line x1="${ml}" y1="${y}" x2="${ml+pw}" y2="${y}" stroke="var(--pos)" stroke-width="1" stroke-dasharray="6,3"/><text x="${ml+pw}" y="${y-4}" font-size="10" fill="var(--pos)" text-anchor="end">dN/dS = 1</text>`; } }

// ── Sliding-window dN/dS along the alignment — positional "Manhattan" ──
if (DATA.window) {
  const box = addSection("dN/dS along the alignment (sliding window)", "winplot");
  const w = DATA.window.filter(p=>p[1]!=null && isFinite(p[1]));
  const W=900,H=360,ml=60,mr=30,mt=20,mb=54,pw=W-ml-mr,ph=H-mt-mb;
  const xmax=Math.max(1,...w.map(p=>p[0])), ymax=Math.max(1.1,...w.map(p=>p[1]))*1.1;
  const X=p=>ml+p/xmax*pw, Y=v=>mt+ph*(1-v/ymax);
  let s=`<svg viewBox="0 0 ${W} ${H}">`;
  for(let i=0;i<=5;i++){ const val=ymax*i/5,y=Y(val);
    s+=`<line x1="${ml}" y1="${y}" x2="${ml+pw}" y2="${y}" stroke="var(--border)" stroke-width="0.5"/><text x="${ml-8}" y="${y}" font-size="10" fill="var(--muted)" text-anchor="end" dominant-baseline="middle">${val.toFixed(2)}</text>`; }
  if(1<=ymax){ const y=Y(1); s+=`<line x1="${ml}" y1="${y}" x2="${ml+pw}" y2="${y}" stroke="var(--pos)" stroke-width="1" stroke-dasharray="6,3"/><text x="${ml+pw}" y="${y-4}" font-size="10" fill="var(--pos)" text-anchor="end">dN/dS = 1</text>`; }
  // connecting line + points
  s+=`<polyline fill="none" stroke="var(--accent)" stroke-width="1.5" opacity="0.6" points="${w.map(p=>X(p[0]).toFixed(1)+','+Y(p[1]).toFixed(1)).join(' ')}"/>`;
  w.forEach(p=>{ const col=p[1]>1?"var(--pos)":"var(--accent)";
    s+=`<circle cx="${X(p[0]).toFixed(1)}" cy="${Y(p[1]).toFixed(1)}" r="3.4" fill="${col}" data-tip="codon ~${p[0]} · dN/dS ${fmt(p[1],3)}"/>`; });
  s+=`<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt+ph}" stroke="var(--fg)"/><line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="var(--fg)"/>`;
  s+=`<text x="${ml+pw/2}" y="${H-6}" font-size="12" fill="var(--muted)" text-anchor="middle">codon position</text>`;
  s+=`<text x="16" y="${mt+ph/2}" font-size="12" fill="var(--muted)" text-anchor="middle" transform="rotate(-90,16,${mt+ph/2})">mean dN/dS</text></svg>`;
  box.innerHTML=s; wireTip(box);
}

// ── Lineage strip-scatter (points per genome + per-lineage mean) ──
if (DATA.lineage) {
  const box = addSection("dN/dS by lineage — points per genome, bar = mean", "linplot");
  const pts = DATA.lineage.filter(d=>d.ratio!=null && isFinite(d.ratio));
  const lins = [...new Set(pts.map(d=>d.lineage))];
  const meanOf={}, nOf={};
  lins.forEach(L=>{ const vs=pts.filter(d=>d.lineage===L).map(d=>d.ratio); nOf[L]=vs.length; meanOf[L]=vs.reduce((a,b)=>a+b,0)/vs.length; });
  const W=Math.max(560, 96*lins.length+140), H=460, ml=60,mr=30,mt=20,mb=110, pw=W-ml-mr, ph=H-mt-mb;
  const ymax=Math.max(1.1, ...pts.map(d=>d.ratio))*1.1;
  const X=i=>ml+(i+0.5)/lins.length*pw, Y=v=>mt+ph*(1-v/ymax);
  let seed=97; const rnd=()=>{seed=(seed*1103515245+12345)&0x7fffffff; return seed/0x7fffffff;};
  // Prefer the mycolorsTB lineage color by name (e.g. L1..L10, A1..A4); else a series slot.
  const lkey=L=>{ const t=String(L).toUpperCase().match(/^(A|L)\s*([0-9]{1,2})/); return t?t[1]+t[2]:null; };
  const lc=L=>{ const k=lkey(L); if(k&&MYCOLORS[k]) return MYCOLORS[k]; const i=lins.indexOf(L); return i<8?`var(${SER[i]})`:"var(--ns)"; };
  const s={v:`<svg viewBox="0 0 ${W} ${H}">`};
  grid(s,ml,pw,mt,ph,ymax,Y,2); neutral(s,ml,pw,ymax,Y);
  lins.forEach((L,i)=>{ const cx=X(i), cw=pw/lins.length*0.62;
    pts.filter(d=>d.lineage===L).forEach(d=>{ const jx=cx+(rnd()-0.5)*cw, cy=Y(d.ratio);
      s.v+=`<circle cx="${jx.toFixed(1)}" cy="${cy.toFixed(1)}" r="3.3" fill="${lc(L)}" opacity="0.6" data-tip="<b>${hesc(d.genome)}</b><br>${hesc(L)} · dN/dS ${fmt(d.ratio,3)}"/>`; });
    const my=Y(meanOf[L]);
    s.v+=`<line x1="${(cx-cw/1.5).toFixed(1)}" y1="${my.toFixed(1)}" x2="${(cx+cw/1.5).toFixed(1)}" y2="${my.toFixed(1)}" stroke="var(--fg)" stroke-width="2.6" data-tip="<b>${hesc(L)}</b><br>mean ${fmt(meanOf[L],3)} (n=${nOf[L]})"/>`;
    s.v+=`<text x="${cx}" y="${mt+ph+14}" font-size="10" fill="var(--muted)" text-anchor="end" transform="rotate(-40,${cx},${mt+ph+14})">${hesc(L)}</text>`; });
  s.v+=`<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt+ph}" stroke="var(--fg)"/><line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="var(--fg)"/>`;
  s.v+=`<text x="16" y="${mt+ph/2}" font-size="12" fill="var(--muted)" text-anchor="middle" transform="rotate(-90,16,${mt+ph/2})">dN/dS</text></svg>`;
  box.innerHTML=s.v+`<div class="legend">`+lins.slice(0,8).map((L,i)=>`<span><i style="background:var(${SER[i]})"></i>${hesc(L)}</span>`).join("")+(lins.length>8?`<span><i style="background:var(--ns)"></i>other</span>`:"")+`</div>`;
  wireTip(box);
}

// ── Group mean ± 95% CI scatter ──
if (DATA.group) {
  const box = addSection("dN/dS by group — mean ± 95% CI", "grpplot");
  const g = DATA.group.filter(d=>d.mean!=null && isFinite(d.mean));
  const W=Math.max(560, 90*g.length+140), H=440, ml=60,mr=30,mt=20,mb=120, pw=W-ml-mr, ph=H-mt-mb;
  const ymax=Math.max(1.1, ...g.map(d=>isFinite(d.ciHigh)?d.ciHigh:d.mean))*1.1;
  const X=i=>ml+(i+0.5)/g.length*pw, Y=v=>mt+ph*(1-v/ymax);
  const s={v:`<svg viewBox="0 0 ${W} ${H}">`};
  grid(s,ml,pw,mt,ph,ymax,Y,2); neutral(s,ml,pw,ymax,Y);
  g.forEach((d,i)=>{ const cx=X(i);
    if(isFinite(d.ciLow)&&isFinite(d.ciHigh)){ s.v+=`<line x1="${cx}" y1="${Y(d.ciLow).toFixed(1)}" x2="${cx}" y2="${Y(d.ciHigh).toFixed(1)}" stroke="var(--muted)" stroke-width="1.5"/>`;
      s.v+=`<line x1="${cx-5}" y1="${Y(d.ciHigh).toFixed(1)}" x2="${cx+5}" y2="${Y(d.ciHigh).toFixed(1)}" stroke="var(--muted)"/><line x1="${cx-5}" y1="${Y(d.ciLow).toFixed(1)}" x2="${cx+5}" y2="${Y(d.ciLow).toFixed(1)}" stroke="var(--muted)"/>`; }
    s.v+=`<circle cx="${cx}" cy="${Y(d.mean).toFixed(1)}" r="4.5" fill="var(--accent)" data-tip="<b>${hesc(d.label)}</b><br>mean ${fmt(d.mean,3)}<br>95% CI [${fmt(d.ciLow,3)}, ${fmt(d.ciHigh,3)}]"/>`;
    s.v+=`<text x="${cx}" y="${mt+ph+14}" font-size="10" fill="var(--muted)" text-anchor="end" transform="rotate(-40,${cx},${mt+ph+14})">${hesc(d.label)}</text>`; });
  s.v+=`<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt+ph}" stroke="var(--fg)"/><line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="var(--fg)"/>`;
  s.v+=`<text x="16" y="${mt+ph/2}" font-size="12" fill="var(--muted)" text-anchor="middle" transform="rotate(-90,16,${mt+ph/2})">mean dN/dS</text></svg>`;
  box.innerHTML=s.v; wireTip(box);
}

// ── dN vs dS scatter (one point per pair) ──
if (DATA.dnds) {
  const box = addSection("dN vs dS — one point per pair (above the line = dN>dS)", "dndsplot");
  const pts = DATA.dnds.filter(p=>p[0]!=null && p[1]!=null && isFinite(p[0]) && isFinite(p[1]));
  const W=560,H=460,ml=64,mr=24,mt=20,mb=54,pw=W-ml-mr,ph=H-mt-mb;
  const amax=Math.max(0.1, ...pts.map(p=>Math.max(p[0],p[1])))*1.08;
  const X=v=>ml+v/amax*pw, Y=v=>mt+ph*(1-v/amax);
  let s=`<svg viewBox="0 0 ${W} ${H}">`;
  for(let i=0;i<=5;i++){ const val=amax*i/5;
    s+=`<line x1="${ml}" y1="${Y(val)}" x2="${ml+pw}" y2="${Y(val)}" stroke="var(--border)" stroke-width="0.5"/><text x="${ml-8}" y="${Y(val)}" font-size="10" fill="var(--muted)" text-anchor="end" dominant-baseline="middle">${val.toFixed(2)}</text>`;
    s+=`<text x="${X(val)}" y="${mt+ph+16}" font-size="10" fill="var(--muted)" text-anchor="middle">${val.toFixed(2)}</text>`; }
  // neutral diagonal dN = dS
  s+=`<line x1="${X(0)}" y1="${Y(0)}" x2="${X(amax)}" y2="${Y(amax)}" stroke="var(--pos)" stroke-width="1" stroke-dasharray="6,3"/><text x="${X(amax)-4}" y="${Y(amax)+14}" font-size="10" fill="var(--pos)" text-anchor="end">dN = dS</text>`;
  pts.forEach(p=>{ const col=p[0]>p[1]?"var(--pos)":"var(--accent)";
    s+=`<circle cx="${X(p[1]).toFixed(1)}" cy="${Y(p[0]).toFixed(1)}" r="2.6" fill="${col}" opacity="0.45" data-tip="dN ${fmt(p[0],3)} · dS ${fmt(p[1],3)} · dN/dS ${p[1]>0?fmt(p[0]/p[1],3):'∞'}"/>`; });
  s+=`<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt+ph}" stroke="var(--fg)"/><line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="var(--fg)"/>`;
  s+=`<text x="${ml+pw/2}" y="${H-6}" font-size="12" fill="var(--muted)" text-anchor="middle">dS</text>`;
  s+=`<text x="16" y="${mt+ph/2}" font-size="12" fill="var(--muted)" text-anchor="middle" transform="rotate(-90,16,${mt+ph/2})">dN</text></svg>`;
  box.innerHTML=s; wireTip(box);
}

// ── Pairwise dN/dS distribution ──
if (DATA.hist) {
  const box = addSection("Pairwise dN/dS distribution", "histplot");
  const h=DATA.hist, cmax=Math.max(1,...h.map(d=>d.count));
  const W=760,H=320,ml=60,mr=20,mt=20,mb=70,pw=W-ml-mr,ph=H-mt-mb,bw=pw/h.length;
  let s=`<svg viewBox="0 0 ${W} ${H}">`;
  for(let i=0;i<=4;i++){ const val=cmax*i/4, y=mt+ph*(1-i/4);
    s+=`<line x1="${ml}" y1="${y}" x2="${ml+pw}" y2="${y}" stroke="var(--border)" stroke-width="0.5"/><text x="${ml-8}" y="${y}" font-size="10" fill="var(--muted)" text-anchor="end" dominant-baseline="middle">${Math.round(val)}</text>`; }
  h.forEach((d,i)=>{ const x=ml+i*bw+3, bh=ph*d.count/cmax, y=mt+ph-bh;
    s+=`<rect x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${(bw-6).toFixed(1)}" height="${bh.toFixed(1)}" fill="var(--accent)" opacity="0.8" data-tip="<b>${hesc(d.label)}</b><br>${d.count} pairs"/>`;
    s+=`<text x="${(x+(bw-6)/2).toFixed(1)}" y="${mt+ph+16}" font-size="9" fill="var(--muted)" text-anchor="middle">${d.label}</text>`; });
  s+=`<line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="var(--fg)"/>`;
  s+=`<text x="${ml+pw/2}" y="${H-6}" font-size="12" fill="var(--muted)" text-anchor="middle">dN/dS bin</text></svg>`;
  box.innerHTML=s; wireTip(box);
}

if(!DATA.lineage && !DATA.group && !DATA.hist && !DATA.dnds && !DATA.window){ sec.innerHTML='<p style="color:var(--muted)">No visualizations available for this run.</p>'; }

// ── Theme toggle + JSON export ────────────────────────────────────────
$("#themeTog").addEventListener("click", ()=>{ const r=document.documentElement;
  const cur=r.getAttribute("data-theme")||(matchMedia("(prefers-color-scheme:dark)").matches?"dark":"light");
  r.setAttribute("data-theme", cur==="dark"?"light":"dark"); });
$("#expJson").addEventListener("click", ()=>{ const b=new Blob([JSON.stringify(DATA,null,1)],{type:"application/json"});
  const a=document.createElement("a"); a.href=URL.createObjectURL(b); a.download="eskaks_dnds.json"; a.click(); URL.revokeObjectURL(a.href); });
"##;
