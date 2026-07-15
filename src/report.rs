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
}

/// JSON-escape a string for embedding in the report.
fn esc(s: &str) -> String {
    let mut o = String::with_capacity(s.len() + 2);
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            '\r' => o.push_str("\\r"),
            '\t' => o.push_str("\\t"),
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
pub fn write_html_report(
    results: &[GenePnPs],
    gw: Option<&GenomeWidePnPs>,
    meta: &ReportMeta,
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

    // ── Embedded data ──────────────────────────────────────────────────────
    let mut data = String::with_capacity(4096 + results.len() * 200);
    data.push_str("{\n");

    // meta
    let ci = meta
        .gw_ci
        .map(|(lo, hi)| format!("[{}, {}]", num(lo), num(hi)))
        .unwrap_or_else(|| "null".to_string());
    let ratio_name = if meta.af_weighted { "πN/πS" } else { "pN/pS" };
    let _ = writeln!(
        data,
        "\"meta\":{{\"samples\":{},\"code\":\"{}\",\"kappa\":{},\"afWeighted\":{},\"fdr\":{},\"minSnps\":{},\"mk\":{},\"mkFixedAf\":{},\"ratioName\":\"{}\"}},",
        meta.n_samples, esc(meta.genetic_code), num(meta.kappa), meta.af_weighted,
        num(meta.fdr), meta.min_snps, meta.mk, num(meta.mk_fixed_af), ratio_name
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
    let _ = writeln!(
        data,
        "\"summary\":{{\"totalGenes\":{},\"genesWithSnps\":{},\"tested\":{},\"significant\":{},\"gwPn\":{},\"gwPs\":{},\"gwRatio\":{},\"gwLabel\":\"{}\",\"gwCi\":{}}},",
        total_genes, genes_with_snps, n_tested, n_sig, gw_pn, gw_ps, gw_ratio, esc(&gw_label), ci
    );

    // genes
    data.push_str("\"genes\":[\n");
    for (i, r) in results.iter().enumerate() {
        let sites = r.n_sites + r.s_sites;
        let exp_n = if sites > 0.0 { r.n_sites / sites } else { f64::NAN };
        let comma = if i + 1 < results.len() { "," } else { "" };
        let _ = write!(
            data,
            "{{\"name\":\"{}\",\"chrom\":\"{}\",\"start\":{},\"end\":{},\"strand\":\"{}\",\"nSites\":{},\"sSites\":{},\"expN\":{},\"pn\":{},\"ps\":{},\"ratio\":{},\"nonsyn\":{},\"syn\":{},\"total\":{},\"p\":{},\"q\":{},\"bonf\":{}",
            esc(&r.name), esc(&r.chrom), r.genome_start, r.genome_end, r.strand,
            num(r.n_sites), num(r.s_sites), num(exp_n), num(r.pn), num(r.ps), num(r.pn_ps),
            num(r.nonsyn_snps), num(r.syn_snps), num(r.total_snps),
            num(r.p_value), num(r.q_value), num(r.p_bonferroni)
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
<title>eskaks pN/pS report</title>
<style>
:root{--bg:#ffffff;--fg:#1a1a1a;--muted:#6b7280;--card:#f6f7f9;--border:#e2e5ea;
--accent:#4a90d9;--pos:#d94a4a;--sig:#d94a4a;--ns:#9aa0a6;--line:#2a7f3f;}
@media (prefers-color-scheme:dark){:root{--bg:#14161a;--fg:#e6e8eb;
--muted:#9aa0a6;--card:#1d2026;--border:#2a2e36;}}
*{box-sizing:border-box}
body{margin:0;font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;
background:var(--bg);color:var(--fg);line-height:1.5}
.wrap{max-width:1100px;margin:0 auto;padding:24px}
h1{font-size:1.5rem;margin:0 0 4px}
.sub{color:var(--muted);font-size:.85rem;margin-bottom:20px}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;margin-bottom:24px}
.card{background:var(--card);border:1px solid var(--border);border-radius:10px;padding:14px}
.card .k{color:var(--muted);font-size:.75rem;text-transform:uppercase;letter-spacing:.04em}
.card .v{font-size:1.35rem;font-weight:650;margin-top:4px}
.card .v.small{font-size:1rem}
section{margin-bottom:28px}
h2{font-size:1.05rem;border-bottom:1px solid var(--border);padding-bottom:6px}
.toolbar{display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin:12px 0}
button.tog{background:var(--card);border:1px solid var(--border);color:var(--fg);
padding:6px 12px;border-radius:8px;cursor:pointer;font-size:.85rem}
button.tog.on{background:var(--accent);color:#fff;border-color:var(--accent)}
input[type=search]{flex:1;min-width:180px;padding:7px 10px;border:1px solid var(--border);
border-radius:8px;background:var(--bg);color:var(--fg);font-size:.9rem}
#plot{width:100%;overflow-x:auto}
svg{max-width:100%;height:auto;display:block}
.tip{position:fixed;pointer-events:none;background:var(--fg);color:var(--bg);
padding:6px 9px;border-radius:6px;font-size:.78rem;opacity:0;transition:opacity .08s;z-index:10;white-space:nowrap}
.tablewrap{overflow-x:auto;border:1px solid var(--border);border-radius:10px}
table{border-collapse:collapse;width:100%;font-size:.83rem}
th,td{padding:7px 10px;text-align:right;white-space:nowrap}
th:first-child,td:first-child{text-align:left;position:sticky;left:0;background:var(--bg)}
thead th{position:sticky;top:0;background:var(--card);cursor:pointer;user-select:none;border-bottom:1px solid var(--border)}
thead th:hover{color:var(--accent)}
thead th.sorted::after{content:" ▲";font-size:.7em}
thead th.sorted.desc::after{content:" ▼"}
tbody tr:nth-child(even){background:var(--card)}
tbody tr.sig td{color:var(--sig);font-weight:600}
.count{color:var(--muted);font-size:.8rem;margin-left:auto}
footer{color:var(--muted);font-size:.75rem;margin-top:24px;text-align:center}
</style>
</head>
<body>
"#;

const BODY: &str = r#"<div class="wrap">
<h1>eskaks — pN/pS report</h1>
<div class="sub" id="meta"></div>
<div class="cards" id="cards"></div>

<section>
<h2>Manhattan</h2>
<div class="toolbar">
  <button class="tog on" id="btnP" data-metric="neglogp">−log10(p)</button>
  <button class="tog" id="btnR" data-metric="ratio">pN/pS</button>
  <span class="count" id="plotNote"></span>
</div>
<div id="plot"></div>
</section>

<section>
<h2>Per-gene table</h2>
<div class="toolbar">
  <input type="search" id="filter" placeholder="Filter by gene / chromosome…">
  <span class="count" id="tableCount"></span>
</div>
<div class="tablewrap"><table id="tbl"><thead></thead><tbody></tbody></table></div>
</section>

<footer>Generated by eskaks · self-contained, no external assets</footer>
</div>
<div class="tip" id="tip"></div>
"#;

const SCRIPT: &str = r##"
const $ = s => document.querySelector(s);
const genes = DATA.genes, S = DATA.summary, M = DATA.meta;
const fmt = (v,d=4) => v==null||!isFinite(v) ? "NA" : Number(v).toFixed(d);
const fmtP = v => v==null||!isFinite(v) ? "NA" : (v!==0 && Math.abs(v)<1e-3 ? v.toExponential(2) : v.toFixed(4));

// ── Meta + summary cards ──────────────────────────────────────────────
$("#meta").textContent =
  `${M.samples} sample(s) · genetic code ${M.code} · kappa ${fmt(M.kappa,2)}` +
  (M.afWeighted ? " · AF-weighted (πN/πS)" : "") +
  (M.minSnps>0 ? ` · min-snps ${M.minSnps}` : "");

function card(k,v,small){return `<div class="card"><div class="k">${k}</div><div class="v${small?" small":""}">${v}</div></div>`;}
let ciTxt = S.gwCi ? ` [${fmt(S.gwCi[0])}, ${fmt(S.gwCi[1])}]` : "";
$("#cards").innerHTML =
  card("Genes analyzed", S.totalGenes) +
  card("Genes with SNPs", S.genesWithSnps) +
  card(`Genome-wide ${M.ratioName}`, fmt(S.gwRatio,4)+`<span style="font-size:.6em;color:var(--muted)">${ciTxt}</span>`, true) +
  card("Selection", S.gwLabel.split(" (")[0], true) +
  card("Genes tested", S.tested) +
  card(`Significant (FDR<${M.fdr})`, S.significant);

// ── Benjamini-Hochberg significance threshold p* ──────────────────────
const ps = genes.map(g=>g.p).filter(v=>v!=null&&isFinite(v)).sort((a,b)=>a-b);
const m = ps.length; let pStar = null;
for (let i=0;i<m;i++){ if (ps[i] <= ((i+1)/m)*M.fdr) pStar = ps[i]; }
const isSig = g => g.q!=null && isFinite(g.q) && g.q < M.fdr;

// ── Interactive Manhattan ─────────────────────────────────────────────
const tip = $("#tip");
let metric = "neglogp";
const W=900,H=420,ml=60,mr=30,mt=20,mb=50,pw=W-ml-mr,ph=H-mt-mb;

function plottable(){ return genes.filter(g => metric==="ratio"
  ? (g.ratio!=null && isFinite(g.ratio) && g.total>0)
  : (g.p!=null && isFinite(g.p))); }

function drawPlot(){
  const pts = plottable();
  const note = $("#plotNote");
  if(!pts.length){ $("#plot").innerHTML=""; note.textContent="no data for this metric"; return; }
  note.textContent = `${pts.length} genes`;
  const yval = g => metric==="ratio" ? g.ratio : -Math.log10(Math.max(g.p,1e-300));
  const xs = pts.map(g=>g.start);
  const xmin=Math.min(...xs), xmax=Math.max(...xs), xr=Math.max(1,xmax-xmin);
  let ymax = Math.max(...pts.map(yval));
  const thrLine = metric==="ratio" ? 1 : (pStar!=null ? -Math.log10(pStar) : null);
  if(thrLine!=null) ymax=Math.max(ymax,thrLine);
  ymax = (ymax||1)*1.1;
  const X = p => ml + ((p-xmin)/xr)*pw;
  const Y = v => mt + ph*(1 - v/ymax);
  let s = `<svg viewBox="0 0 ${W} ${H}" role="img">`;
  s += `<rect x="0" y="0" width="${W}" height="${H}" fill="none"/>`;
  for(let i=0;i<=5;i++){ const val=ymax*i/5, y=Y(val);
    s+=`<line x1="${ml}" y1="${y}" x2="${ml+pw}" y2="${y}" stroke="var(--border)" stroke-width="0.5"/>`;
    s+=`<text x="${ml-8}" y="${y}" font-size="10" fill="var(--muted)" text-anchor="end" dominant-baseline="middle">${val.toFixed(metric==="ratio"?2:1)}</text>`; }
  if(thrLine!=null){ const y=Y(thrLine);
    s+=`<line x1="${ml}" y1="${y}" x2="${ml+pw}" y2="${y}" stroke="var(--line)" stroke-width="1" stroke-dasharray="6,3"/>`;
    s+=`<text x="${ml+pw}" y="${y-4}" font-size="10" fill="var(--line)" text-anchor="end">${metric==="ratio"?"pN/pS = 1":"BH FDR "+M.fdr}</text>`; }
  pts.forEach((g,i)=>{ const r=Math.min(8,Math.max(2.5,Math.sqrt(g.total||1)));
    const color = isSig(g) ? "var(--sig)" : (metric==="ratio" ? (g.ratio<1?"var(--accent)":"var(--pos)") : "var(--ns)");
    s+=`<circle data-i="${genes.indexOf(g)}" cx="${X(g.start).toFixed(1)}" cy="${Y(yval(g)).toFixed(1)}" r="${r.toFixed(1)}" fill="${color}" opacity="0.78"/>`; });
  s+=`<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt+ph}" stroke="var(--fg)" stroke-width="1.2"/>`;
  s+=`<line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="var(--fg)" stroke-width="1.2"/>`;
  s+=`<text x="${ml+pw/2}" y="${H-10}" font-size="12" fill="var(--muted)" text-anchor="middle">Genome position</text>`;
  s+=`<text x="16" y="${mt+ph/2}" font-size="12" fill="var(--muted)" text-anchor="middle" transform="rotate(-90,16,${mt+ph/2})">${metric==="ratio"?"pN/pS":"−log10(p)"}</text>`;
  s+=`</svg>`;
  $("#plot").innerHTML = s;
}

$("#plot").addEventListener("mousemove", e=>{
  const t=e.target;
  if(t.tagName==="circle"){ const g=genes[+t.dataset.i];
    tip.innerHTML = `<b>${g.name}</b> (${g.chrom}:${g.start} ${g.strand})<br>pN/pS ${fmt(g.ratio,3)} · p ${fmtP(g.p)} · q ${fmtP(g.q)}<br>${g.nonsyn}N / ${g.syn}S SNPs`;
    tip.style.opacity=1; tip.style.left=(e.clientX+12)+"px"; tip.style.top=(e.clientY+12)+"px";
  } else { tip.style.opacity=0; }
});
$("#plot").addEventListener("mouseleave", ()=>tip.style.opacity=0);
document.querySelectorAll("button.tog").forEach(b=>b.addEventListener("click",()=>{
  metric=b.dataset.metric;
  document.querySelectorAll("button.tog").forEach(x=>x.classList.toggle("on",x===b));
  drawPlot();
}));

// ── Sortable / filterable table ───────────────────────────────────────
let cols = [
  ["name","Gene"],["chrom","Chrom"],["start","Start"],["strand","±"],
  ["nSites","N sites"],["sSites","S sites"],["ratio",M.ratioName],
  ["nonsyn","N SNPs"],["syn","S SNPs"],["expN","Exp N frac"],
  ["p","p"],["q","q (BH)"],["bonf","p (Bonf)"]
];
if(M.mk) cols=cols.concat([["dn","Dn"],["ds","Ds"],["pnMk","Pn"],["psMk","Ps"],["ni","NI"],["alpha","α"],["fisherP","Fisher p"]]);
const pcols = new Set(["p","q","bonf","fisherP"]);
let sortKey="p", sortDesc=false;

$("#tbl thead").innerHTML = "<tr>"+cols.map(c=>`<th data-k="${c[0]}">${c[1]}</th>`).join("")+"</tr>";
function renderTable(){
  const f = $("#filter").value.toLowerCase();
  let rows = genes.filter(g => !f || g.name.toLowerCase().includes(f) || (g.chrom||"").toLowerCase().includes(f));
  rows.sort((a,b)=>{ let x=a[sortKey], y=b[sortKey];
    const xn=(x==null||(typeof x==="number"&&!isFinite(x))), yn=(y==null||(typeof y==="number"&&!isFinite(y)));
    if(xn&&yn) return 0; if(xn) return 1; if(yn) return -1;
    if(typeof x==="string") return sortDesc ? y.localeCompare(x) : x.localeCompare(y);
    return sortDesc ? y-x : x-y; });
  $("#tbl tbody").innerHTML = rows.map(g=>"<tr"+(isSig(g)?' class="sig"':"")+">"+cols.map(c=>{
    let v=g[c[0]];
    let txt = typeof v==="string" ? v : (c[0]==="start"||["nonsyn","syn","total","dn","ds","pnMk","psMk"].includes(c[0]) ? (v==null?"NA":v) : (pcols.has(c[0])?fmtP(v):fmt(v, c[0]==="ratio"||c[0]==="ni"||c[0]==="alpha"?4:2)));
    return `<td>${txt}</td>`;
  }).join("")+"</tr>").join("");
  document.querySelectorAll("#tbl thead th").forEach(th=>{
    th.classList.toggle("sorted", th.dataset.k===sortKey);
    th.classList.toggle("desc", th.dataset.k===sortKey && sortDesc); });
  $("#tableCount").textContent = `${rows.length} / ${genes.length} genes`;
}
$("#tbl thead").addEventListener("click", e=>{ const th=e.target.closest("th"); if(!th)return;
  const k=th.dataset.k; if(k===sortKey) sortDesc=!sortDesc; else {sortKey=k; sortDesc=false;} renderTable(); });
$("#filter").addEventListener("input", renderTable);

drawPlot(); renderTable();
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
<div class="cards" id="cards"></div>
<div id="sections"></div>
<footer>Generated by eskaks · self-contained, no external assets</footer>
</div>
<div class="tip" id="tip"></div>
"#;

const FASTA_SCRIPT: &str = r##"
const $ = s => document.querySelector(s);
const M = DATA.meta, tip = $("#tip");
const fmt = (v,d=4) => v==null||!isFinite(v) ? "NA" : Number(v).toFixed(d);

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
  const s={v:`<svg viewBox="0 0 ${W} ${H}">`};
  grid(s,ml,pw,mt,ph,ymax,Y,2); neutral(s,ml,pw,ymax,Y);
  lins.forEach((L,i)=>{ const cx=X(i), cw=pw/lins.length*0.62;
    pts.filter(d=>d.lineage===L).forEach(d=>{ const jx=cx+(rnd()-0.5)*cw, cy=Y(d.ratio);
      s.v+=`<circle cx="${jx.toFixed(1)}" cy="${cy.toFixed(1)}" r="3.3" fill="var(--accent)" opacity="0.55" data-tip="<b>${d.genome}</b><br>${L} · dN/dS ${fmt(d.ratio,3)}"/>`; });
    const my=Y(meanOf[L]);
    s.v+=`<line x1="${(cx-cw/1.5).toFixed(1)}" y1="${my.toFixed(1)}" x2="${(cx+cw/1.5).toFixed(1)}" y2="${my.toFixed(1)}" stroke="var(--fg)" stroke-width="2.6" data-tip="<b>${L}</b><br>mean ${fmt(meanOf[L],3)} (n=${nOf[L]})"/>`;
    s.v+=`<text x="${cx}" y="${mt+ph+14}" font-size="10" fill="var(--muted)" text-anchor="end" transform="rotate(-40,${cx},${mt+ph+14})">${L}</text>`; });
  s.v+=`<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt+ph}" stroke="var(--fg)"/><line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="var(--fg)"/>`;
  s.v+=`<text x="16" y="${mt+ph/2}" font-size="12" fill="var(--muted)" text-anchor="middle" transform="rotate(-90,16,${mt+ph/2})">dN/dS</text></svg>`;
  box.innerHTML=s.v; wireTip(box);
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
    s.v+=`<circle cx="${cx}" cy="${Y(d.mean).toFixed(1)}" r="4.5" fill="var(--accent)" data-tip="<b>${d.label}</b><br>mean ${fmt(d.mean,3)}<br>95% CI [${fmt(d.ciLow,3)}, ${fmt(d.ciHigh,3)}]"/>`;
    s.v+=`<text x="${cx}" y="${mt+ph+14}" font-size="10" fill="var(--muted)" text-anchor="end" transform="rotate(-40,${cx},${mt+ph+14})">${d.label}</text>`; });
  s.v+=`<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt+ph}" stroke="var(--fg)"/><line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="var(--fg)"/>`;
  s.v+=`<text x="16" y="${mt+ph/2}" font-size="12" fill="var(--muted)" text-anchor="middle" transform="rotate(-90,16,${mt+ph/2})">mean dN/dS</text></svg>`;
  box.innerHTML=s.v; wireTip(box);
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
    s+=`<rect x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${(bw-6).toFixed(1)}" height="${bh.toFixed(1)}" fill="var(--accent)" opacity="0.8" data-tip="<b>${d.label}</b><br>${d.count} pairs"/>`;
    s+=`<text x="${(x+(bw-6)/2).toFixed(1)}" y="${mt+ph+16}" font-size="9" fill="var(--muted)" text-anchor="middle">${d.label}</text>`; });
  s+=`<line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="var(--fg)"/>`;
  s+=`<text x="${ml+pw/2}" y="${H-6}" font-size="12" fill="var(--muted)" text-anchor="middle">dN/dS bin</text></svg>`;
  box.innerHTML=s; wireTip(box);
}

if(!DATA.lineage && !DATA.group && !DATA.hist){ sec.innerHTML='<p style="color:var(--muted)">Run with --lineage, --group-average, or the default pairwise mode to populate visualizations.</p>'; }
"##;
