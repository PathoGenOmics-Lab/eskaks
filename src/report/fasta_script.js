
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
