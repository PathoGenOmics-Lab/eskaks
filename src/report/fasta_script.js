
const $ = s => document.querySelector(s);
const M = DATA.meta, tip = $("#tip");
// HTML-escape genome/lineage/group names before putting them in SVG text, tooltips
// or legends, so a name with < > & " can't break markup or inject HTML.
const hesc = s => String(s==null?"":s).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;");
const fmt = (v,d=4) => v==null||!isFinite(v) ? "NA" : Number(v).toFixed(d);
// Tip HTML must NOT round-trip through a data-* attribute: the HTML parser decodes the
// attribute value once, undoing hesc() and letting a name like `<img onerror=…>` run
// when the tip is later assigned via innerHTML. Stash the HTML in an array and put only
// its index in the attribute; wireTip reads it back unchanged.
const TIPS=[];
const tipData = h => `data-tip="${TIPS.push(h)-1}"`;
// Max over an array WITHOUT spreading it into call arguments (a large alignment has
// O(n^2) pairs / many windows, and the spread overflows the call-argument limit).
const maxOf = (arr, seed, f) => { let m=seed; for(const x of arr){ const v=f?f(x):x; if(v>m)m=v; } return m; };
// PathoGenOmics-Lab mycolorsTB — M. tuberculosis lineage colors (keyed by name).
const MYCOLORS = {A1:"#d1ae00",A2:"#8ef5c8",A3:"#73c2ff",A4:"#ff9cdb",L1:"#ff3091",L2:"#001aff",
  L3:"#8a0bd2",L4:"#ff0000",L5:"#995200",L6:"#1eb040",L7:"#fbff00",L8:"#ff9d00",L9:"#37ff30",L10:"#8fbda1"};
const SER = ["--s1","--s2","--s3","--s4","--s5","--s6","--s7","--s8"];

// ── "i" interpretation help — same popover mechanism as the vcf report ──────────
const HELP = {
  pairs:{title:"Pairs & valid pairs",lead:"How many pairwise sequence comparisons were formed, and how many yielded a usable dN/dS.",x:"",y:"",
    read:["Pairs: every sequence pair compared (n·(n−1)/2 for n aligned sequences).","Valid pairs: those with dS &gt; 0 and a finite dN/dS; a pair with no synonymous differences is undefined and excluded."],
    watch:["A wide gap between pairs and valid pairs means many comparisons had zero synonymous substitutions — common for near-identical or very short sequences."]},
  pooled:{title:"Pooled dN/dS",lead:"One overall ratio pooled across all valid pairs: Σ dN / Σ dS, so more-divergent pairs weigh more.",x:"",y:"",
    read:["Below 1: purifying selection — nonsynonymous change is suppressed, the usual signal for conserved coding sequence.","Above 1: excess nonsynonymous divergence — candidate positive/diversifying selection.","This is a point estimate over pairs; the report reports no per-ratio significance test or CI."],
    watch:["Pooling is dominated by the most divergent pairs; it summarises the set, not any single lineage.","dN/dS between very closely related sequences is unstable — read it alongside the per-pair scatter and the sliding window."]},
  meandnds:{title:"Mean dN / dS",lead:"The mean of per-pair dN and of per-pair dS, reported separately (not the mean of the ratios).",x:"",y:"",
    read:["Compare with the pooled ratio: if the two diverge, a few high-divergence pairs are driving the pooled value."],
    watch:["Averaging dN and dS separately avoids the instability of averaging ratios when some pairs have tiny dS."]},
  window:{title:"dN/dS along the alignment (sliding window)",lead:"Mean dN/dS in a window slid across codon positions, revealing where along the sequence selective pressure varies.",x:"Codon position along the alignment (window centre).",y:"Mean dN/dS in the window; the dashed line marks the neutral expectation dN/dS = 1.",
    read:["Dips well below 1: strongly constrained regions (active sites, structural cores).","Peaks above 1 (red): a local excess of amino-acid change — candidate diversifying sites.","The connecting line traces the positional trend; isolated points are single windows."],
    watch:["Short windows are noisy; a lone peak from few substitutions is fragile.","The window track shares the model and alignment of the pooled estimate — a misalignment shifts the whole curve."]},
  lineage:{title:"dN/dS by lineage",lead:"Per-genome dN/dS grouped by lineage, with a bar at each lineage's mean, showing between-lineage variation in selective pressure.",x:"Lineage (grouped along x).",y:"dN/dS for each genome (points); the bar is the lineage mean. Dashed line = neutral (1).",
    read:["A lineage sitting well below 1: pervasive purifying selection.","A lineage shifted upward: relaxed constraint or diversifying selection worth a closer look.","Spread within a lineage shows how consistent the signal is across its genomes."],
    watch:["Lineages with few genomes give unstable means.","Lineage colours follow the mycolorsTB scheme by name (L1–L10, A1–A4) where recognised, else a series colour."]},
  group:{title:"dN/dS by group — mean ± 95% CI",lead:"Group mean dN/dS with a 95% confidence interval, so a group whose interval excludes 1 shows a supported departure from neutrality.",x:"Group.",y:"Mean dN/dS; whiskers = 95% CI. Dashed line = neutral (1).",
    read:["Interval entirely below 1: well-supported purifying selection.","Interval spanning 1: not distinguishable from neutral at this sample size.","Interval entirely above 1: supported diversifying selection."],
    watch:["The CI reflects between-pair scatter and group size; small groups give wide intervals.","Overlapping intervals between groups mean the difference is not resolved."]},
  dndsscatter:{title:"dN vs dS — one point per pair",lead:"Each pair plotted by synonymous (dS) and nonsynonymous (dN) divergence; the diagonal is neutrality (dN = dS).",x:"dS — synonymous substitutions per synonymous site.",y:"dN — nonsynonymous substitutions per nonsynonymous site.",
    read:["Below the diagonal (dN &lt; dS): purifying selection — the majority signal for coding sequence.","Above the diagonal (dN &gt; dS, red): excess amino-acid change — candidate positive selection.","Distance out along the diagonal reflects overall divergence between the pair."],
    watch:["Points near the origin are near-identical pairs where dN/dS is very unstable.","A pair at dS ≈ 0 has an undefined ratio and is excluded from the pooled estimate."]},
  histdist:{title:"Pairwise dN/dS distribution",lead:"Histogram of per-pair dN/dS, showing the dominant selective regime across the pair set and any tail of high-ratio pairs.",x:"Per-pair dN/dS bin (1.0 = neutral).",y:"Number of pairs in each bin.",
    read:["Mass below 1: purifying selection dominates the set.","A bin around 1: pairs indistinguishable from neutral, often low-divergence.","A right tail beyond 1: candidate diversifying pairs — inspect them in the dN-vs-dS scatter."],
    watch:["Low-divergence pairs give unstable ratios that smear the distribution; a bin's position is not a significance test."]},
};
function helpHtml(h){
  return '<h4>'+h.title+'</h4><p class="lead">'+h.lead+'</p>'
    +(h.x?'<div class="ax"><b>x</b> — '+h.x+'</div>':'')
    +(h.y?'<div class="ax"><b>y</b> — '+h.y+'</div>':'')
    +((h.read&&h.read.length)?'<ul>'+h.read.map(function(s){return '<li>'+s+'</li>';}).join('')+'</ul>':'')
    +((h.watch&&h.watch.length)?'<div class="warn"><b>Watch out</b> — '+h.watch.join(' · ')+'</div>':'');
}
const infoBtn = key => (HELP[key]?`<button class="info" data-help="${key}" title="What is this?" aria-label="What is this?">i</button>`:"");
const infopop=$("#infopop"); let openInfo=null;
function positionPop(btn){
  const r=btn.getBoundingClientRect(), pw=infopop.offsetWidth, ph=infopop.offsetHeight, m=8;
  const left=Math.max(m, Math.min(r.left, window.innerWidth-pw-m));
  let top=r.bottom+6; if(top+ph>window.innerHeight-m) top=r.top-ph-6;   // flip above if it would overflow below
  top=Math.max(m, Math.min(top, window.innerHeight-ph-m));               // clamp fully into the viewport
  infopop.style.left=left+'px'; infopop.style.top=top+'px';
}
function showInfo(btn){ const h=HELP[btn.dataset.help]; if(!h||!infopop) return;
  infopop.innerHTML=helpHtml(h); infopop.hidden=false; positionPop(btn);
  btn.classList.add('on'); openInfo=btn; }
function hideInfo(){ if(!infopop) return; infopop.hidden=true; if(openInfo)openInfo.classList.remove('on'); openInfo=null; }
document.addEventListener('click', function(e){ const b=e.target.closest&&e.target.closest('.info');
  if(b){ e.stopPropagation(); if(openInfo===b){hideInfo();}else{hideInfo();showInfo(b);} return; }
  if(!(e.target.closest&&e.target.closest('#infopop'))) hideInfo(); });
document.addEventListener('keydown', function(e){ if(e.key==='Escape' && openInfo){ hideInfo(); e.stopImmediatePropagation(); } });
window.addEventListener('resize', hideInfo);
window.addEventListener('scroll', hideInfo, true);

$("#meta").textContent = `model ${M.model} · ${M.totalPairs} pairs (${M.validPairs} valid)`;
const card = (k,v,cls,help) => `<div class="card${cls?" "+cls:""}"><div class="k">${k}${help?infoBtn(help):""}</div><div class="v">${v}</div></div>`;
$("#cards").innerHTML =
  card("Pairs", M.totalPairs, "", "pairs") + card("Valid pairs", M.validPairs, "", "pairs") +
  card("Pooled dN/dS", fmt(M.pooled), "hero", "pooled") + card("Mean dN / dS", fmt(M.meanDn,3)+" / "+fmt(M.meanDs,3), "", "meandnds");

// ── Hero verdict banner: plain-language direction of the pooled ratio (mirrors the vcf report).
// This is a pooled point estimate over lineage/group pairs, with no per-ratio significance test,
// so the wording states the direction only and never claims significance.
(function renderVerdict(){
  const v=$("#verdict"); if(!v) return;
  const gw=M.pooled; let lab, dotc;
  if(gw==null||!isFinite(gw)){ lab="dN/dS undetermined"; dotc="var(--ns)"; }
  else if(gw<1){ lab="purifying selection (dN/dS < 1)"; dotc="var(--accent)"; }
  else if(gw>1){ lab="positive / diversifying selection (dN/dS > 1)"; dotc="var(--pos)"; }
  else { lab="neutral (dN/dS = 1)"; dotc="var(--ns)"; }
  const val=(gw==null||!isFinite(gw))?"":fmt(gw,3)+" ";
  v.innerHTML=
    '<span class="dot" style="background:'+dotc+';color:'+dotc+'"></span>'+
    '<div style="flex:1;min-width:0"><div class="vlab">Pooled dN/dS '+val+'— '+lab+'</div>'+
    '<div class="vsub">pooled across '+M.validPairs+' of '+M.totalPairs+' pair(s) · '+hesc(M.model)+' · point estimate</div></div>'+
    infoBtn("pooled");
  v.style.setProperty("--vc", dotc);   // regime colour drives the banner's accent edge + wash
})();

const sec = $("#sections");
function addSection(title, id, help){ const d=document.createElement("section");
  d.innerHTML=`<h2><span class="htitle">${title}</span>${help?infoBtn(help):""}</h2><div style="overflow-x:auto" id="${id}"></div>`;
  sec.appendChild(d); return $("#"+id); }
function wireTip(c){
  c.addEventListener("mousemove", e=>{ const info=e.target.getAttribute && e.target.getAttribute("data-tip");
    if(info){ tip.innerHTML=TIPS[+info]||""; tip.style.opacity=1; tip.style.left=(e.clientX+12)+"px"; tip.style.top=(e.clientY+12)+"px"; } else tip.style.opacity=0; });
  c.addEventListener("mouseleave", ()=>tip.style.opacity=0);
}
function grid(s,ml,pw,mt,ph,ymax,Y,dec){ for(let i=0;i<=5;i++){ const val=ymax*i/5,y=Y(val);
  s.v+=`<line x1="${ml}" y1="${y}" x2="${ml+pw}" y2="${y}" stroke="var(--border)" stroke-width="0.5"/><text x="${ml-8}" y="${y}" font-size="10" fill="var(--muted)" text-anchor="end" dominant-baseline="middle">${val.toFixed(dec)}</text>`; } }
function neutral(s,ml,pw,ymax,Y){ if(1<=ymax){ const y=Y(1);
  s.v+=`<line x1="${ml}" y1="${y}" x2="${ml+pw}" y2="${y}" stroke="var(--pos)" stroke-width="1" stroke-dasharray="6,3"/><text x="${ml+pw}" y="${y-4}" font-size="10" fill="var(--pos)" text-anchor="end">dN/dS = 1</text>`; } }

// ── Sliding-window dN/dS along the alignment — positional "Manhattan" ──
if (DATA.window) {
  const box = addSection("dN/dS along the alignment (sliding window)", "winplot", "window");
  const w = DATA.window.filter(p=>p[1]!=null && isFinite(p[1]));
  const W=900,H=360,ml=60,mr=30,mt=20,mb=54,pw=W-ml-mr,ph=H-mt-mb;
  const xmax=maxOf(w,1,p=>p[0]), ymax=maxOf(w,1.1,p=>p[1])*1.1;
  const X=p=>ml+p/xmax*pw, Y=v=>mt+ph*(1-v/ymax);
  let s=`<svg viewBox="0 0 ${W} ${H}">`;
  for(let i=0;i<=5;i++){ const val=ymax*i/5,y=Y(val);
    s+=`<line x1="${ml}" y1="${y}" x2="${ml+pw}" y2="${y}" stroke="var(--border)" stroke-width="0.5"/><text x="${ml-8}" y="${y}" font-size="10" fill="var(--muted)" text-anchor="end" dominant-baseline="middle">${val.toFixed(2)}</text>`; }
  if(1<=ymax){ const y=Y(1); s+=`<line x1="${ml}" y1="${y}" x2="${ml+pw}" y2="${y}" stroke="var(--pos)" stroke-width="1" stroke-dasharray="6,3"/><text x="${ml+pw}" y="${y-4}" font-size="10" fill="var(--pos)" text-anchor="end">dN/dS = 1</text>`; }
  // connecting line + points
  s+=`<polyline fill="none" stroke="var(--accent)" stroke-width="1.5" opacity="0.6" points="${w.map(p=>X(p[0]).toFixed(1)+','+Y(p[1]).toFixed(1)).join(' ')}"/>`;
  w.forEach(p=>{ const col=p[1]>1?"var(--pos)":"var(--accent)";
    s+=`<circle cx="${X(p[0]).toFixed(1)}" cy="${Y(p[1]).toFixed(1)}" r="3.4" fill="${col}" ${tipData(`codon ~${p[0]} · dN/dS ${fmt(p[1],3)}`)}/>`; });
  s+=`<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt+ph}" stroke="var(--fg)"/><line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="var(--fg)"/>`;
  s+=`<text x="${ml+pw/2}" y="${H-6}" font-size="12" fill="var(--muted)" text-anchor="middle">codon position</text>`;
  s+=`<text x="16" y="${mt+ph/2}" font-size="12" fill="var(--muted)" text-anchor="middle" transform="rotate(-90,16,${mt+ph/2})">mean dN/dS</text></svg>`;
  box.innerHTML=s; wireTip(box);
}

// ── Lineage strip-scatter (points per genome + per-lineage mean) ──
if (DATA.lineage) {
  const box = addSection("dN/dS by lineage — points per genome, bar = mean", "linplot", "lineage");
  const pts = DATA.lineage.filter(d=>d.ratio!=null && isFinite(d.ratio));
  const lins = [...new Set(pts.map(d=>d.lineage))];
  const meanOf={}, nOf={};
  lins.forEach(L=>{ const vs=pts.filter(d=>d.lineage===L).map(d=>d.ratio); nOf[L]=vs.length; meanOf[L]=vs.reduce((a,b)=>a+b,0)/vs.length; });
  const W=Math.max(560, 96*lins.length+140), H=460, ml=60,mr=30,mt=20,mb=110, pw=W-ml-mr, ph=H-mt-mb;
  const ymax=maxOf(pts,1.1,d=>d.ratio)*1.1;
  const X=i=>ml+(i+0.5)/lins.length*pw, Y=v=>mt+ph*(1-v/ymax);
  let seed=97; const rnd=()=>{seed=(seed*1103515245+12345)&0x7fffffff; return seed/0x7fffffff;};
  // Prefer the mycolorsTB lineage color by name (e.g. L1..L10, A1..A4); else a series slot.
  const lkey=L=>{ const t=String(L).toUpperCase().match(/^(A|L)\s*([0-9]{1,2})/); return t?t[1]+t[2]:null; };
  const lc=L=>{ const k=lkey(L); if(k&&MYCOLORS[k]) return MYCOLORS[k]; const i=lins.indexOf(L); return i<8?`var(${SER[i]})`:"var(--ns)"; };
  const s={v:`<svg viewBox="0 0 ${W} ${H}">`};
  grid(s,ml,pw,mt,ph,ymax,Y,2); neutral(s,ml,pw,ymax,Y);
  lins.forEach((L,i)=>{ const cx=X(i), cw=pw/lins.length*0.62;
    pts.filter(d=>d.lineage===L).forEach(d=>{ const jx=cx+(rnd()-0.5)*cw, cy=Y(d.ratio);
      s.v+=`<circle cx="${jx.toFixed(1)}" cy="${cy.toFixed(1)}" r="3.3" fill="${lc(L)}" opacity="0.6" ${tipData(`<b>${hesc(d.genome)}</b><br>${hesc(L)} · dN/dS ${fmt(d.ratio,3)}`)}/>`; });
    const my=Y(meanOf[L]);
    s.v+=`<line x1="${(cx-cw/1.5).toFixed(1)}" y1="${my.toFixed(1)}" x2="${(cx+cw/1.5).toFixed(1)}" y2="${my.toFixed(1)}" stroke="var(--fg)" stroke-width="2.6" ${tipData(`<b>${hesc(L)}</b><br>mean ${fmt(meanOf[L],3)} (n=${nOf[L]})`)}/>`;
    s.v+=`<text x="${cx}" y="${mt+ph+14}" font-size="10" fill="var(--muted)" text-anchor="end" transform="rotate(-40,${cx},${mt+ph+14})">${hesc(L)}</text>`; });
  s.v+=`<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt+ph}" stroke="var(--fg)"/><line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="var(--fg)"/>`;
  s.v+=`<text x="16" y="${mt+ph/2}" font-size="12" fill="var(--muted)" text-anchor="middle" transform="rotate(-90,16,${mt+ph/2})">dN/dS</text></svg>`;
  box.innerHTML=s.v+`<div class="legend">`+lins.slice(0,8).map(L=>`<span><i style="background:${lc(L)}"></i>${hesc(L)}</span>`).join("")+(lins.length>8?`<span><i style="background:var(--ns)"></i>other</span>`:"")+`</div>`;
  wireTip(box);
}

// ── Group mean ± 95% CI scatter ──
if (DATA.group) {
  const box = addSection("dN/dS by group — mean ± 95% CI", "grpplot", "group");
  const g = DATA.group.filter(d=>d.mean!=null && isFinite(d.mean));
  const W=Math.max(560, 90*g.length+140), H=440, ml=60,mr=30,mt=20,mb=120, pw=W-ml-mr, ph=H-mt-mb;
  const ymax=maxOf(g,1.1,d=>isFinite(d.ciHigh)?d.ciHigh:d.mean)*1.1;
  const X=i=>ml+(i+0.5)/g.length*pw, Y=v=>mt+ph*(1-v/ymax);
  const s={v:`<svg viewBox="0 0 ${W} ${H}">`};
  grid(s,ml,pw,mt,ph,ymax,Y,2); neutral(s,ml,pw,ymax,Y);
  g.forEach((d,i)=>{ const cx=X(i);
    if(isFinite(d.ciLow)&&isFinite(d.ciHigh)){ s.v+=`<line x1="${cx}" y1="${Y(d.ciLow).toFixed(1)}" x2="${cx}" y2="${Y(d.ciHigh).toFixed(1)}" stroke="var(--muted)" stroke-width="1.5"/>`;
      s.v+=`<line x1="${cx-5}" y1="${Y(d.ciHigh).toFixed(1)}" x2="${cx+5}" y2="${Y(d.ciHigh).toFixed(1)}" stroke="var(--muted)"/><line x1="${cx-5}" y1="${Y(d.ciLow).toFixed(1)}" x2="${cx+5}" y2="${Y(d.ciLow).toFixed(1)}" stroke="var(--muted)"/>`; }
    s.v+=`<circle cx="${cx}" cy="${Y(d.mean).toFixed(1)}" r="4.5" fill="var(--accent)" ${tipData(`<b>${hesc(d.label)}</b><br>mean ${fmt(d.mean,3)}<br>95% CI [${fmt(d.ciLow,3)}, ${fmt(d.ciHigh,3)}]`)}/>`;
    s.v+=`<text x="${cx}" y="${mt+ph+14}" font-size="10" fill="var(--muted)" text-anchor="end" transform="rotate(-40,${cx},${mt+ph+14})">${hesc(d.label)}</text>`; });
  s.v+=`<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt+ph}" stroke="var(--fg)"/><line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="var(--fg)"/>`;
  s.v+=`<text x="16" y="${mt+ph/2}" font-size="12" fill="var(--muted)" text-anchor="middle" transform="rotate(-90,16,${mt+ph/2})">mean dN/dS</text></svg>`;
  box.innerHTML=s.v; wireTip(box);
}

// ── dN vs dS scatter (one point per pair) ──
if (DATA.dnds) {
  const box = addSection("dN vs dS — one point per pair (above the line = dN>dS)", "dndsplot", "dndsscatter");
  const pts = DATA.dnds.filter(p=>p[0]!=null && p[1]!=null && isFinite(p[0]) && isFinite(p[1]));
  const W=560,H=460,ml=64,mr=24,mt=20,mb=54,pw=W-ml-mr,ph=H-mt-mb;
  const amax=maxOf(pts,0.1,p=>Math.max(p[0],p[1]))*1.08;
  const X=v=>ml+v/amax*pw, Y=v=>mt+ph*(1-v/amax);
  let s=`<svg viewBox="0 0 ${W} ${H}">`;
  for(let i=0;i<=5;i++){ const val=amax*i/5;
    s+=`<line x1="${ml}" y1="${Y(val)}" x2="${ml+pw}" y2="${Y(val)}" stroke="var(--border)" stroke-width="0.5"/><text x="${ml-8}" y="${Y(val)}" font-size="10" fill="var(--muted)" text-anchor="end" dominant-baseline="middle">${val.toFixed(2)}</text>`;
    s+=`<text x="${X(val)}" y="${mt+ph+16}" font-size="10" fill="var(--muted)" text-anchor="middle">${val.toFixed(2)}</text>`; }
  // neutral diagonal dN = dS
  s+=`<line x1="${X(0)}" y1="${Y(0)}" x2="${X(amax)}" y2="${Y(amax)}" stroke="var(--pos)" stroke-width="1" stroke-dasharray="6,3"/><text x="${X(amax)-4}" y="${Y(amax)+14}" font-size="10" fill="var(--pos)" text-anchor="end">dN = dS</text>`;
  pts.forEach(p=>{ const col=p[0]>p[1]?"var(--pos)":"var(--accent)";
    s+=`<circle cx="${X(p[1]).toFixed(1)}" cy="${Y(p[0]).toFixed(1)}" r="2.6" fill="${col}" opacity="0.45" ${tipData(`dN ${fmt(p[0],3)} · dS ${fmt(p[1],3)} · dN/dS ${p[1]>0?fmt(p[0]/p[1],3):'∞'}`)}/>`; });
  s+=`<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt+ph}" stroke="var(--fg)"/><line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="var(--fg)"/>`;
  s+=`<text x="${ml+pw/2}" y="${H-6}" font-size="12" fill="var(--muted)" text-anchor="middle">dS</text>`;
  s+=`<text x="16" y="${mt+ph/2}" font-size="12" fill="var(--muted)" text-anchor="middle" transform="rotate(-90,16,${mt+ph/2})">dN</text></svg>`;
  box.innerHTML=s; wireTip(box);
}

// ── Pairwise dN/dS distribution ──
if (DATA.hist) {
  const box = addSection("Pairwise dN/dS distribution", "histplot", "histdist");
  const h=DATA.hist, cmax=maxOf(h,1,d=>d.count);
  const W=760,H=320,ml=60,mr=20,mt=20,mb=70,pw=W-ml-mr,ph=H-mt-mb,bw=pw/h.length;
  let s=`<svg viewBox="0 0 ${W} ${H}">`;
  for(let i=0;i<=4;i++){ const val=cmax*i/4, y=mt+ph*(1-i/4);
    s+=`<line x1="${ml}" y1="${y}" x2="${ml+pw}" y2="${y}" stroke="var(--border)" stroke-width="0.5"/><text x="${ml-8}" y="${y}" font-size="10" fill="var(--muted)" text-anchor="end" dominant-baseline="middle">${Math.round(val)}</text>`; }
  h.forEach((d,i)=>{ const x=ml+i*bw+3, bh=ph*d.count/cmax, y=mt+ph-bh;
    s+=`<rect x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${(bw-6).toFixed(1)}" height="${bh.toFixed(1)}" fill="var(--accent)" opacity="0.8" ${tipData(`<b>${hesc(d.label)}</b><br>${d.count} pairs`)}/>`;
    s+=`<text x="${(x+(bw-6)/2).toFixed(1)}" y="${mt+ph+16}" font-size="9" fill="var(--muted)" text-anchor="middle">${hesc(d.label)}</text>`; });
  s+=`<line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="var(--fg)"/>`;
  s+=`<text x="${ml+pw/2}" y="${H-6}" font-size="12" fill="var(--muted)" text-anchor="middle">dN/dS bin</text></svg>`;
  box.innerHTML=s; wireTip(box);
}

if(!DATA.lineage && !DATA.group && !DATA.hist && !DATA.dnds && !DATA.window){ sec.innerHTML='<p style="color:var(--muted)">No visualizations available for this run.</p>'; }

// ── Theme toggle · CSV / JSON export · print ──────────────────────────
$("#themeTog").addEventListener("click", ()=>{ const r=document.documentElement;
  const cur=r.getAttribute("data-theme")||(matchMedia("(prefers-color-scheme:dark)").matches?"dark":"light");
  r.setAttribute("data-theme", cur==="dark"?"light":"dark"); });

function dl(name,txt,type){ const b=new Blob([txt],{type}); const a=document.createElement("a");
  a.href=URL.createObjectURL(b); a.download=name; a.click(); URL.revokeObjectURL(a.href); }
// Numbers pass through verbatim (full precision, and never a formula); only STRING cells are
// quoted and, if they open with =+-@ (tab/CR), prefixed with ' to defuse spreadsheet formula injection.
const csvCell = v => {
  if(typeof v==="number") return isFinite(v)?String(v):"";
  let s=String(v==null?"":v);
  // Defuse spreadsheet formula injection: =/@ always start a formula; +/- only when more than a
  // lone sign follows (so a bare "+"/"-" datum is left intact, not turned into "'+").
  if(/^[=@\t\r\n]/.test(s) || (s.length>1 && /^[+\-]/.test(s))) s="'"+s;
  return /[",\n\r]/.test(s) ? '"'+s.replace(/"/g,'""')+'"' : s;
};
const ratio = (dn,ds) => (dn!=null && isFinite(dn) && ds>0) ? dn/ds : "";   // undefined unless dN is estimable and dS>0
// Export the richest labelled table this run carries (lineage > group > per-pair),
// mirroring how the JSON export dumps whatever DATA holds.
function buildCsv(){
  let rows;
  if(DATA.lineage && DATA.lineage.length){
    rows=[["genome","lineage","dN/dS"]];
    DATA.lineage.forEach(d=>rows.push([d.genome, d.lineage, d.ratio]));
  } else if(DATA.group && DATA.group.length){
    rows=[["group","mean_dN/dS","ci_low","ci_high"]];
    DATA.group.forEach(d=>rows.push([d.label, d.mean, d.ciLow, d.ciHigh]));
  } else if(DATA.dnds && DATA.dnds.length){
    rows=[["pair","dN","dS","dN/dS"]];
    DATA.dnds.forEach((p,i)=>rows.push([i+1, p[0], p[1], ratio(p[0],p[1])]));
  } else {
    rows=[["metric","value"],["model",M.model],["pairs",M.totalPairs],["valid_pairs",M.validPairs],
      ["pooled_dN/dS",M.pooled],["mean_dN",M.meanDn],["mean_dS",M.meanDs]];
  }
  return rows.map(r=>r.map(csvCell).join(",")).join("\n");
}
$("#expCsv").addEventListener("click", ()=>dl("eskaks_dnds.csv", buildCsv(), "text/csv"));
$("#expJson").addEventListener("click", ()=>dl("eskaks_dnds.json", JSON.stringify(DATA,null,1), "application/json"));
$("#printBtn").addEventListener("click", ()=>{ hideInfo(); window.print(); });
