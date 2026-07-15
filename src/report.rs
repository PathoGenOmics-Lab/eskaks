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
