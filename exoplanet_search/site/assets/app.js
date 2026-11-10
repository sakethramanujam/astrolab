const C={acc:'#7cc4ff',acc2:'#ffb86b',ok:'#5eead4',warn:'#fbbf24',bad:'#f87171',sub:'#c084fc',muted:'#8b93b3',grid:'rgba(124,196,255,.08)'};
const fmt=(x,d=0)=>x==null?'–':Number(x).toLocaleString(undefined,{maximumFractionDigits:d,minimumFractionDigits:d});
const el=(t,a={},h='')=>{const e=document.createElement(t);Object.assign(e,a);if(h)e.innerHTML=h;return e};
const rel=iso=>{if(!iso)return'';const s=(Date.now()-Date.parse(iso))/1000;
  if(s<90)return'just now';if(s<5400)return Math.round(s/60)+' min ago';if(s<172800)return Math.round(s/3600)+' h ago';return Math.round(s/86400)+' days ago'};
if(window.Chart){Chart.defaults.color='#c3c9e3';Chart.defaults.font.family='Inter,system-ui,sans-serif';Chart.defaults.borderColor=C.grid;
  Chart.defaults.plugins.legend.labels.boxWidth=10;}

fetch('data/progress.json',{cache:'no-store'}).then(r=>r.json()).then(render).catch(e=>{
  document.getElementById('updated').textContent='Could not load data/progress.json: '+e});

function render(d){
  const t=d.totals;
  document.getElementById('updated').textContent=
    `last updated ${new Date(d.generated_at).toUTCString()} (${rel(d.generated_at)}) · pipeline v${d.pipeline_version}`+
    (t.planet_archive_crossmatch?'':' · NASA Exoplanet Archive cross-match pending (archive unavailable, no cached copy yet)')+
    Object.entries(d.catalogs||{}).map(([k,v])=>v.fetched_at?` · ${k.replace('.csv','')} ${v.fetched_at.slice(0,10)}`:'').join('');
  const rate=t.toi_hosts?100*t.toi_recovered/t.toi_hosts:0;
  const stats=[
    ['Stars searched',fmt(t.stars_searched),`${fmt(t.errors)} errors`,'acc'],
    ['Sectors',`${t.sectors_done}<span style="font-size:18px;color:var(--muted)"> / ${t.sectors_available}</span>`,`${t.sectors_with_results} with results · next: ${(d.queue_next||[]).join(', ')||'–'}`,''],
    ['Known TOIs re-found',`${fmt(rate,0)}%`,`${t.toi_recovered} of ${t.toi_hosts} TOI hosts`,'ok'],
    ['Needs review',fmt(t.needs_review),'passed every screen, no catalogue match','acc2'],
    ['Submitted',fmt(t.submitted),'to Planet Hunters TESS / ExoFOP','sub']];
  const S=document.getElementById('stats');
  for(const [k,v,s,c] of stats)S.append(el('div',{className:'card stat'},`<div class="k">${k}</div><div class="v ${c}">${v}</div><div class="s">${s}</div>`));

  // queue
  const total=t.sectors_available||1, done=t.sectors_done, prog=d.sectors.filter(s=>s.status==='in_progress').length;
  document.getElementById('qbar').innerHTML=`<i style="width:${100*done/total}%;background:${C.ok}"></i><i style="width:${Math.max(100*prog/total,prog?1.2:0)}%;background:${C.acc2}"></i>`;
  document.getElementById('qlegend').innerHTML=`<span><span class="dot" style="background:${C.ok}"></span>${done} done</span><span><span class="dot" style="background:${C.acc2}"></span>${prog} in progress</span><span><span class="dot" style="background:rgba(124,196,255,.25)"></span>${total-done-prog} queued</span>`;
  const st=document.getElementById('sectors');
  st.innerHTML='<thead><tr><th>Sector</th><th>Status</th><th class="num">Searched</th><th class="num">Passed</th><th class="num">TOIs re-found</th><th class="num">Review</th></tr></thead>';
  const tb=el('tbody');
  for(const s of [...d.sectors].sort((a,b)=>b.sector-a.sector)){
    const pct=s.n_targets?` <span class="kbd">/ ${fmt(s.n_targets)}</span>`:'';
    tb.append(el('tr',{},`<td>S${s.sector}</td><td>${s.status==='done'?'✅ done':'⏳ in progress'}</td><td class="num">${fmt(s.n_searched)}${pct}</td><td class="num">${s.n_passed}</td><td class="num">${s.toi_recovered} / ${s.toi_hosts}</td><td class="num">${s.needs_review}</td>`));
  }
  st.append(tb);

  if(!window.Chart)return renderTables(d);
  new Chart(cRecovery,{type:'doughnut',data:{labels:['re-found & passed','found, failed a screen','missed'],
    datasets:[{data:[t.toi_recovered,t.toi_period_found-t.toi_recovered,t.toi_hosts-t.toi_period_found],backgroundColor:[C.ok,C.warn,'rgba(139,147,179,.35)'],borderWidth:0}]},
    options:{maintainAspectRatio:false,cutout:'68%',plugins:{legend:{position:'right'}}}});
  const fu=d.deep_funnel;
  if(fu){
    new Chart(cFunnel,{type:'bar',data:{labels:['passed screens, no match','deep-vetted','on-target (no blend)','not binary-like','survive everything'],
      datasets:[{data:[fu.unmatched_passed,fu.deep_vetted,fu.survive_blend,fu.survive_eb,fu.survive_all],
      backgroundColor:['rgba(124,196,255,.45)','rgba(124,196,255,.65)',C.acc,C.sub,C.acc2],borderWidth:0}]},
      options:{indexAxis:'y',maintainAspectRatio:false,plugins:{legend:{display:false}}}});
    const dn={centroid_offset:'dip centroid off target',dip_off_target:'darkest pixel off target',unresolved_neighbour:'unresolved bright neighbour',
      neighbour_could_mimic:'neighbour could mimic (centroid inconclusive)',centroid_inconclusive:'centroid inconclusive',centroid_shift_minor:'small centroid shift (warning)',known_eb:'known eclipsing binary',
      v_shape:'V-shaped (grazing binary-like)',secondary_scan:'secondary eclipse found',rp_gt_2rjup:'too big (> 2 Jupiter radii)'};
    const de=Object.entries(fu.by_flag).filter(x=>x[1]>0).sort((a,b)=>b[1]-a[1]);
    new Chart(cDeep,{type:'bar',data:{labels:de.map(x=>dn[x[0]]||x[0]),datasets:[{data:de.map(x=>x[1]),backgroundColor:'rgba(192,132,252,.6)',borderWidth:0}]},
      options:{indexAxis:'y',maintainAspectRatio:false,plugins:{legend:{display:false}}}});
  }
  const sh=d.charts.sde_hist, ph=d.charts.period_hist;
  new Chart(cSde,{type:'bar',data:{labels:sh.edges.slice(0,-1).map(x=>x.toFixed(0)),datasets:[{label:'stars',data:sh.counts,
    backgroundColor:sh.edges.slice(0,-1).map(x=>x>=9?C.acc2:'rgba(124,196,255,.55)'),borderWidth:0,barPercentage:1,categoryPercentage:1}]},
    options:{maintainAspectRatio:false,plugins:{legend:{display:false}},scales:{y:{type:'logarithmic',title:{display:true,text:'stars (log)'}},x:{title:{display:true,text:'SDE (orange = above threshold 9)'}}}}});
  new Chart(cPeriod,{type:'bar',data:{labels:ph.edges.slice(0,-1).map(x=>x<2?x.toFixed(2):x.toFixed(1)),datasets:[{label:'stars',data:ph.counts,backgroundColor:'rgba(192,132,252,.6)',borderWidth:0,barPercentage:1,categoryPercentage:1}]},
    options:{maintainAspectRatio:false,plugins:{legend:{display:false}},scales:{x:{title:{display:true,text:'period (days, log bins) — pile-ups are spacecraft systematics'}}}}});
  const fl=Object.entries(d.charts.flags).sort((a,b)=>b[1]-a[1]);
  const names={low_sde:'weak periodogram peak',low_snr:'low transit SNR',too_long:'too long for the star',few_transits:'< 3 transits',too_deep:'too deep (> 5%)',
    odd_even:'odd ≠ even depth',variable:'variable star',secondary:'secondary eclipse',rp_too_large:'too big for a planet',common_period:'shared systematic period',no_dip:'no dip'};
  new Chart(cFlags,{type:'bar',data:{labels:fl.map(f=>names[f[0]]||f[0]),datasets:[{data:fl.map(f=>f[1]),backgroundColor:'rgba(248,113,113,.55)',borderWidth:0}]},
    options:{indexAxis:'y',maintainAspectRatio:false,plugins:{legend:{display:false}}}});
  new Chart(cSky,{type:'scatter',data:{datasets:[{data:d.charts.sky.map(([x,y])=>({x,y})),pointRadius:1.2,backgroundColor:'rgba(124,196,255,.6)'}]},
    options:{maintainAspectRatio:false,plugins:{legend:{display:false}},scales:{x:{min:0,max:360,reverse:true,title:{display:true,text:'RA (deg)'}},y:{min:-90,max:90,title:{display:true,text:'Dec (deg)'}}}}});
  new Chart(cRuns,{type:'bar',data:{labels:d.runs.map(r=>`S${r.sector} · ${r.started_at.slice(0,10)}`),datasets:[
    {label:'stars searched',data:d.runs.map(r=>r.n_searched),backgroundColor:'rgba(124,196,255,.6)',yAxisID:'y'},
    {type:'line',label:'median seconds per star',data:d.runs.map(r=>r.median_star_s),borderColor:C.acc2,backgroundColor:C.acc2,yAxisID:'y1'}]},
    options:{maintainAspectRatio:false,scales:{y:{position:'left'},y1:{position:'right',grid:{display:false}}}}});
  renderTables(d);
}

function renderTables(d){
  const rt=document.getElementById('runs');
  rt.innerHTML='<thead><tr><th>Started (UTC)</th><th>Sector</th><th class="num">Stars</th><th class="num">Errors</th><th class="num">Wall</th><th class="num">s / star</th><th>Backend</th></tr></thead>';
  const rb=el('tbody');
  for(const r of [...d.runs].reverse())rb.append(el('tr',{},`<td>${r.started_at.replace('T',' ').replace('+00:00','')}</td><td>S${r.sector}</td><td class="num">${fmt(r.n_searched)}</td><td class="num">${r.n_errors}</td><td class="num">${fmt(r.search_wall_s/60,1)} min</td><td class="num">${fmt(r.median_star_s,2)}</td><td>${r.backend}${r.workers?` ×${r.workers}`:''}</td>`));
  rt.append(rb);

  const states=['all','needs_review','known','submitted','rejected'], label={all:'All',needs_review:'Needs review',known:'Known',submitted:'Submitted',rejected:'Rejected'};
  let cur='all', key='rank', dir=-1;
  const RANK={submitted:4,needs_review:3,known:2,rejected:1};
  d.candidates.forEach(c=>c.rank=RANK[c.status]*1000+(c.sde||0));
  const chips=document.getElementById('chips');
  for(const s of states){const n=s==='all'?d.candidates.length:d.candidates.filter(c=>c.status===s).length;
    const b=el('button',{className:'chip'+(s===cur?' on':''),textContent:`${label[s]} · ${n}`});
    b.onclick=()=>{cur=s;[...chips.children].forEach(x=>x.classList.remove('on'));b.classList.add('on');draw()};chips.append(b)}
  const cols=[['tic','TIC',0],['sector','Sector',0],['status','Status'],['known','Match'],['period','Period (d)',4],['duration_h','Dur (h)',1],['depth_ppm','Depth (ppm)',0],
    ['sde','SDE',1],['snr','SNR',1],['n_transits','Transits',0],['rp_earth','Rp (R⊕)',1],['tmag','Tmag',1],['allflags','Flags']];
  d.candidates.forEach(c=>c.allflags=[c.flags,c.deep_flags].filter(Boolean).join(',')+(c.deep_vetted?'':(c.status==='needs_review'?' (deep vetting pending)':'')));
  const tab=document.getElementById('cands');
  function draw(){
    const rows=d.candidates.filter(c=>cur==='all'||c.status===cur).sort((a,b)=>{const x=a[key],y=b[key];return (x==null)-(y==null)||(x>y?1:x<y?-1:0)*dir});
    tab.innerHTML='';const th=el('thead');const tr=el('tr');
    for(const [k,l,dg] of cols){const h=el('th',{className:dg!=null?'num':'',textContent:l+(k===key?(dir<0?' ↓':' ↑'):'')});h.onclick=()=>{dir=key===k?-dir:-1;key=k;draw()};tr.append(h)}
    th.append(tr);tab.append(th);const tb=el('tbody');
    for(const c of rows){const r=el('tr',{className:c.plot?'click':''});
      for(const [k,,dg] of cols){let v=c[k];
        if(k==='status')v=`<span class="badge b-${v}">${label[v]||v}</span>`;
        else if(k==='known')v=v||(c.toi_host?`<span class="kbd">TOI ${c.toi_host} host, other period</span>`:'–');
        else if(k==='allflags')v=`<span class="flags">${v||'—'}</span>`;
        else if(k==='tic')v=`<a href="${c.exofop}" target="_blank" rel="noopener" onclick="event.stopPropagation()">${v}</a>`;
        else if(dg!=null)v=fmt(v,dg);
        r.append(el('td',{className:dg!=null&&k!=='tic'?'num':''},v))}
      if(c.plot)r.onclick=()=>openPlot(c);tb.append(r)}
    tab.append(tb)}
  draw();
}
function openPlot(c){
  dlgT.textContent=`TIC ${c.tic} · sector ${c.sector} · ${c.known||'no catalogue match'}`;
  dlgI.src=c.plot;
  dlgB.innerHTML=`Period ${fmt(c.period,5)} d · depth ${fmt(c.depth_ppm)} ppm · duration ${fmt(c.duration_h,1)} h · SDE ${fmt(c.sde,1)} · SNR ${fmt(c.snr,1)} · odd/even ${fmt(c.odd_even_sigma,1)}σ · secondary ${fmt(c.secondary_sigma,1)}σ<br>
  Screens: <span class="flags">${c.flags||'all passed'}</span><br>
  Deep vetting: ${c.deep_vetted?`<span class="flags">${c.deep_flags||'all passed'}</span> · dip centroid offset ${fmt(c.centroid_offset_px,2)} px (${fmt(c.centroid_offset_sigma,1)}σ) · Gaia neighbours that could mimic: ${c.gaia_n_can_mimic??'–'}${c.gaia_mimic_min_sep!=null?` (closest ${fmt(c.gaia_mimic_min_sep,1)}″)`:''} · V-ratio ${fmt(c.v_ratio,2)} · best secondary ${fmt(c.sec_scan_sigma,1)}σ · ${fmt(c.rp_rjup,2)} R<sub>J</sub>${c.eb_match?` · <b>${c.eb_match}</b>`:''}`:'not run yet'}<br> · <a href="${c.exofop}" target="_blank" rel="noopener">ExoFOP page</a> · panels: raw + trend, periodogram, fold, odd/even, phase 0.5.`;
  dlg.showModal()}
