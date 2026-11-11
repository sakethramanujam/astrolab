import { el, fmt, slider, toggle, starfield } from './ui.js';

const C = { acc: '#7cc4ff', acc2: '#ffb86b', ok: '#5eead4', warn: '#fbbf24', bad: '#f87171', sub: '#c084fc', muted: '#8b93b3', grid: 'rgba(124,196,255,.08)' };
const rel = (iso) => {
  if (!iso) return '';
  const s = (Date.now() - Date.parse(iso)) / 1000;
  if (s < 90) return 'just now'; if (s < 5400) return Math.round(s / 60) + ' min ago';
  if (s < 172800) return Math.round(s / 3600) + ' h ago'; return Math.round(s / 86400) + ' days ago';
};
const LABEL = { all: 'All', needs_review: 'Needs review', uncalibrated: 'Not calibrated', known: 'Known', submitted: 'Submitted', rejected: 'Rejected' };
const FLAG_TXT = {
  centroid_offset: 'dip centroid off target', dip_off_target: 'darkest pixel off target', unresolved_neighbour: 'unresolved bright neighbour',
  neighbour_could_mimic: 'neighbour could mimic', centroid_inconclusive: 'centroid inconclusive', centroid_shift_minor: 'small centroid shift (warning)',
  known_eb: 'known eclipsing binary', v_shape: 'V-shaped', secondary_scan: 'secondary eclipse found', rp_gt_2rjup: 'too big (> 2 R_J)',
};
const REJ = new Set(['centroid_offset', 'dip_off_target', 'known_eb', 'v_shape', 'secondary_scan', 'rp_gt_2rjup', 'low_sde', 'low_snr', 'few_transits', 'odd_even', 'secondary', 'too_deep', 'too_long', 'variable', 'rp_too_large', 'common_period']);
starfield();
if (window.Chart) {
  Chart.defaults.color = '#c3c9e3'; Chart.defaults.font.family = 'Inter,system-ui,sans-serif'; Chart.defaults.borderColor = C.grid;
  Chart.defaults.plugins.legend.labels.boxWidth = 10;
  if (matchMedia('(prefers-reduced-motion: reduce)').matches) Chart.defaults.animation = false;
}

fetch('data/progress.json', { cache: 'no-store' }).then((r) => r.json()).then(render).catch((e) => {
  document.getElementById('updated').textContent = 'Could not load data/progress.json: ' + e;
});

function render(d) {
  const t = d.totals;
  document.getElementById('updated').textContent =
    `last updated ${new Date(d.generated_at).toUTCString()} (${rel(d.generated_at)}) · pipeline v${d.pipeline_version}` +
    (t.planet_archive_crossmatch ? '' : ' · NASA Exoplanet Archive cross-match pending (archive unavailable, no cached copy yet)') +
    Object.entries(d.catalogs || {}).map(([k, v]) => (v.fetched_at ? ` · ${k.replace('.csv', '')} ${v.fetched_at.slice(0, 10)}` : '')).join('');
  const nh = t.toi_hosts_calibrated ?? t.toi_hosts;
  const rate = nh ? (100 * t.toi_recovered) / nh : 0;
  const S = document.getElementById('stats');
  for (const [k, v, s, c] of [
    ['Stars searched', fmt(t.stars_searched), `${fmt(t.errors)} errors`, 'acc'],
    ['Sectors', `${t.sectors_done}<span style="font-size:18px;color:var(--muted)"> / ${t.sectors_available}</span>`, `${t.sectors_with_results} with results · next: ${(d.queue_next || []).join(', ') || '–'}`, ''],
    ['Known TOIs re-found', `${fmt(rate, 0)}%`, `${t.toi_recovered} of ${nh} TOI hosts in calibrated (v1) rows, at S ≥ interim S*`, 'ok'],
    ['Needs review', fmt(t.needs_review), `passed every check at S ≥ interim S*, no catalogue match; not detections · ${fmt(t.uncalibrated ?? 0)} more not calibrated`, 'acc2'],
    ['Submitted', fmt(t.submitted), 'to Planet Hunters TESS / ExoFOP', 'sub'],
  ]) S.append(el('div', { class: 'card stat', html: `<div class="k">${k}</div><div class="v ${c}">${v}</div><div class="s">${s}</div>` }));

  const table = candidates(d);
  sectorMap(d, table);
  sectorTable(d);
  charts(d);
  runs(d);
}

/* ───── sector map: one tile per TESS sector ───── */
function sectorMap(d, table) {
  const q = d.queue?.length ? d.queue : d.sectors.map((s) => ({ sector: s.sector, status: s.status, n_targets: s.n_targets, n_searched: s.n_searched }));
  const bySector = Object.fromEntries(d.sectors.map((s) => [s.sector, s]));
  const next = new Set(d.queue_next || []);
  const YEARS = [1, 14, 27, 40, 56, 70, 84, 97, 110, 123, 136];
  const root = document.getElementById('smap');
  const info = document.getElementById('sinfo');
  const cnt = { done: 0, in_progress: 0, pending: 0 };
  q.forEach((s) => (cnt[s.status] = (cnt[s.status] || 0) + 1));
  document.getElementById('qlegend').innerHTML =
    `<span><span class="dot" style="background:${C.ok}"></span>${cnt.done} done</span>` +
    `<span><span class="dot" style="background:${C.acc2}"></span>${cnt.in_progress} in progress</span>` +
    `<span><span class="dot" style="background:rgba(124,196,255,.25)"></span>${cnt.pending} queued</span>` +
    `<span><span class="dot" style="background:${C.warn}"></span>has signals to review</span>`;
  const show = (s) => {
    const r = bySector[s.sector];
    info.innerHTML = [['Sector', `S${s.sector}`], ['Status', s.status.replace('_', ' ') + (next.has(s.sector) && s.status !== 'done' ? ' (up next)' : '')],
      ['Stars searched', `${fmt(s.n_searched)}${s.n_targets ? ' / ' + fmt(s.n_targets) : ''}`],
      ['Passed (S ≥ S*)', r ? (r.calibrated ? fmt(r.n_passed) : 'not calibrated') : '–'], ['TOIs re-found', r ? `${r.toi_recovered} / ${r.toi_hosts}` : '–'], ['Needs review', r ? fmt(r.needs_review) : '–']]
      .map(([k, v]) => `<div><span>${k}</span><b>${v}</b></div>`).join('');
  };
  for (let i = 0; i < YEARS.length - 1; i++) {
    const ss = q.filter((s) => s.sector >= YEARS[i] && s.sector < YEARS[i + 1]);
    if (!ss.length) continue;
    const tiles = el('div', { class: 'tiles' });
    for (const s of ss) {
      const r = bySector[s.sector];
      const frac = s.n_targets ? Math.min(1, s.n_searched / s.n_targets) : 0;
      const label = `Sector ${s.sector}: ${s.status.replace('_', ' ')}${s.n_targets ? `, ${fmt(s.n_searched)} of ${fmt(s.n_targets)} stars` : ''}${r ? `, ${r.needs_review} need review` : ''}`;
      const b = el('button', { type: 'button', class: `tile ${s.status}${next.has(s.sector) && s.status !== 'done' ? ' next' : ''}`, 'aria-label': label, title: label,
        'data-s': s.sector },
        s.status === 'in_progress' ? el('span', { class: 'fill', style: `height:${(frac * 100).toFixed(0)}%` }) : null,
        r && r.needs_review ? el('span', { class: 'rv', 'aria-hidden': 'true' }) : null, el('b', {}, s.sector));
      b.addEventListener('mouseenter', () => show(s));
      b.addEventListener('focus', () => show(s));
      b.addEventListener('click', () => {
        if (!r) { show(s); return; }
        const on = !b.classList.contains('on');
        root.querySelectorAll('.tile.on').forEach((x) => x.classList.remove('on'));
        b.classList.toggle('on', on);
        b.setAttribute('aria-pressed', String(on));
        table.setSector(on ? s.sector : '');
        if (on) document.getElementById('candidates').scrollIntoView({ behavior: matchMedia('(prefers-reduced-motion: reduce)').matches ? 'auto' : 'smooth' });
      });
      tiles.append(b);
    }
    root.append(el('div', { class: 'yr' }, el('span', {}, `Year ${i + 1}`), tiles));
  }
  const latest = q.find((s) => s.status === 'in_progress') || [...q].reverse().find((s) => s.status === 'done') || q[0];
  if (latest) show(latest);
}

const short = (n) => (n == null ? '–' : n >= 1e4 ? `${(n / 1e3).toFixed(1)}k` : fmt(n));
/** Build a responsive table: real table semantics on desktop, label/value cards below 720 px. */
function rtable(table, cols, rows) {
  table.classList.add('rt');
  table.setAttribute('role', 'table');
  const thead = el('thead', { role: 'rowgroup' }, el('tr', { role: 'row' }, cols.map(([, l, cls]) => el('th', { scope: 'col', role: 'columnheader', class: cls || '' }, l))));
  const tbody = el('tbody', { role: 'rowgroup' }, rows.map((r) => el('tr', { role: 'row' }, cols.map(([k, l, cls], i) =>
    el('td', { role: 'cell', class: cls || '', 'data-label': l, html: `<span class="v">${r[i] ?? '–'}</span>` })))));
  table.replaceChildren(...[...table.querySelectorAll('caption')], thead, tbody);
}

function sectorTable(d) {
  rtable(document.getElementById('sectors'),
    [['s', 'Sector'], ['st', 'Status'], ['n', 'Searched', 'num'], ['p', 'Passed', 'num'], ['t', 'TOIs found', 'num'], ['r', 'Review', 'num']],
    [...d.sectors].sort((a, b) => b.sector - a.sector).map((x) => [`S${x.sector}`, x.status === 'done' ? '✅ done' : '⏳ running',
      `${short(x.n_searched)}${x.n_targets && x.n_targets !== x.n_searched ? ` / ${short(x.n_targets)}` : ''}`, fmt(x.n_passed), `${x.toi_recovered} / ${x.toi_hosts}`, fmt(x.needs_review)]));
}

/* ───── candidate explorer: search, filter, sort; state lives in the URL hash ───── */
function candidates(d) {
  const all = d.candidates;
  const RANK = { submitted: 5, needs_review: 4, uncalibrated: 3, known: 2, rejected: 1 };
  all.forEach((c) => {
    c.rank = RANK[c.status] * 1000 + (c.sde || 0);
    c.allflags = [c.flags, c.deep_flags].filter(Boolean).join(',');
    c.search = `${c.tic} ${c.known || ''} ${c.toi_host ? 'toi ' + c.toi_host : ''} s${c.sector} ${c.allflags}`.toLowerCase();
  });
  const flagSet = [...new Set(all.flatMap((c) => c.allflags.split(',').filter(Boolean)))].sort();
  const sectors = [...new Set(all.map((c) => c.sector))].sort((a, b) => b - a);
  const hp = new URLSearchParams(location.hash.slice(1));
  const st = {
    status: hp.get('status') || 'all', q: hp.get('q') || '', sector: hp.get('sector') || '', flag: hp.get('flag') || '',
    pmin: +hp.get('pmin') || 0.3, pmax: +hp.get('pmax') || 30, rmax: +hp.get('rmax') || 40, tmax: +hp.get('tmax') || 17, sdemin: +hp.get('sde') || 0,
    key: hp.get('sort') || 'rank', dir: +hp.get('dir') || -1,
  };
  const save = () => {
    const p = new URLSearchParams();
    if (st.status !== 'all') p.set('status', st.status);
    for (const [k, v, dflt] of [['q', st.q, ''], ['sector', st.sector, ''], ['flag', st.flag, ''], ['pmin', st.pmin, 0.3], ['pmax', st.pmax, 30],
      ['rmax', st.rmax, 40], ['tmax', st.tmax, 17], ['sde', st.sdemin, 0], ['sort', st.key, 'rank'], ['dir', st.dir, -1]])
      if (String(v) !== String(dflt)) p.set(k, typeof v === 'number' ? +v.toFixed(3) : v);
    history.replaceState(null, '', p.toString() ? '#' + p : location.pathname + location.search);
  };

  const F = document.getElementById('filters');
  const search = el('input', { type: 'search', placeholder: 'TIC, TOI, flag…', 'aria-label': 'Search candidates', value: st.q });
  search.addEventListener('input', () => { st.q = search.value.trim().toLowerCase(); draw(); });
  const secSel = el('select', { 'aria-label': 'Sector' }, el('option', { value: '' }, 'All sectors'), sectors.map((s) => el('option', { value: s }, `Sector ${s}`)));
  secSel.value = st.sector;
  secSel.addEventListener('change', () => { st.sector = secSel.value; draw(); });
  const flagSel = el('select', { 'aria-label': 'Flag' }, el('option', { value: '' }, 'Any flags'), el('option', { value: '-' }, 'No flags at all'),
    flagSet.map((f) => el('option', { value: f }, FLAG_TXT[f] || f.replace(/_/g, ' '))));
  flagSel.value = st.flag;
  flagSel.addEventListener('change', () => { st.flag = flagSel.value; draw(); });
  const sP0 = slider({ id: 'pmin', label: 'Period ≥', min: 0.3, max: 30, log: true, value: st.pmin, fmt: (v) => `${v.toFixed(2)} d`, onInput: (v) => { st.pmin = v; draw(); } });
  const sP1 = slider({ id: 'pmax', label: 'Period ≤', min: 0.3, max: 30, log: true, value: st.pmax, fmt: (v) => `${v.toFixed(2)} d`, onInput: (v) => { st.pmax = v; draw(); } });
  const sR = slider({ id: 'rmax', label: 'Radius ≤', min: 0.5, max: 40, log: true, value: st.rmax, fmt: (v) => `${v.toFixed(1)} R⊕`, onInput: (v) => { st.rmax = v; draw(); } });
  const sT = slider({ id: 'tmax', label: 'Tmag ≤ (brighter)', min: 5, max: 17, step: 0.5, value: st.tmax, fmt: (v) => v.toFixed(1), onInput: (v) => { st.tmax = v; draw(); } });
  const sS = slider({ id: 'sde', label: 'SDE ≥', min: 0, max: 30, step: 0.5, value: st.sdemin, fmt: (v) => v.toFixed(1), onInput: (v) => { st.sdemin = v; draw(); } });
  const reset = el('button', { type: 'button', class: 'btn' }, 'Reset filters');
  reset.addEventListener('click', () => {
    Object.assign(st, { status: 'all', q: '', sector: '', flag: '', pmin: 0.3, pmax: 30, rmax: 40, tmax: 17, sdemin: 0 });
    search.value = ''; secSel.value = ''; flagSel.value = '';
    sP0.set(0.3); sP1.set(30); sR.set(40); sT.set(17); sS.set(0); chipsSync();
    document.querySelectorAll('.tile.on').forEach((x) => x.classList.remove('on'));
    draw();
  });
  F.append(el('div', { class: 'ctl' }, el('span', { class: 'lbl' }, 'Search'), search), el('div', { class: 'ctl' }, el('span', { class: 'lbl' }, 'Sector'), secSel),
    el('div', { class: 'ctl' }, el('span', { class: 'lbl' }, 'Flags'), flagSel), sS.el, sP0.el, sP1.el, sR.el, sT.el, el('div', { class: 'ctl' }, reset));

  const chips = document.getElementById('chips');
  const chipBtns = ['all', 'needs_review', 'uncalibrated', 'known', 'submitted', 'rejected'].map((s) => {
    const b = el('button', { type: 'button', class: 'chip', 'data-s': s }, `${LABEL[s]} · ${s === 'all' ? all.length : all.filter((c) => c.status === s).length}`);
    b.addEventListener('click', () => { st.status = s; chipsSync(); draw(); });
    chips.append(b);
    return b;
  });
  const chipsSync = () => chipBtns.forEach((b) => { const on = b.dataset.s === st.status; b.classList.toggle('on', on); b.setAttribute('aria-pressed', String(on)); });
  chipsSync();

  // Key columns only; everything else lives in an expandable detail row (and on the candidate page).
  const pct = (v) => (v == null ? '–' : v >= 10000 ? `${(v / 1e4).toFixed(v >= 1e5 ? 0 : 1)}%` : `${fmt(v)} ppm`);
  const per = (v) => (v == null ? '–' : v >= 10 ? v.toFixed(2) : v.toFixed(3));
  const flagChips = (c) => {
    const f = c.allflags ? c.allflags.split(',') : [];
    if (!f.length) return c.status === 'needs_review' && !c.deep_vetted ? '<span class="fchip">deep vetting pending</span>' : '<span class="kbd">none</span>';
    return f.map((x) => `<span class="fchip${REJ.has(x) ? ' bad' : ''}">${FLAG_TXT[x] || x.replace(/_/g, ' ')}</span>`).join('');
  };
  const cols = [
    ['tic', 'TIC', 'c-tic'], ['status', 'Status', 'c-st'], ['period', 'Period (d)', 'num c-p'], ['depth_ppm', 'Depth', 'num c-d'],
    ['rp_earth', 'Radius (R⊕)', 'num c-r'], ['sde', 'SDE', 'num c-s'], ['allflags', 'Flags', 'c-f'],
  ];
  const tab = document.getElementById('cands');
  tab.classList.add('rt', 'cands');
  tab.setAttribute('role', 'table');
  const thead = el('thead', { role: 'rowgroup' }), tbody = el('tbody', { role: 'rowgroup' });
  const htr = el('tr', { role: 'row' }, el('th', { class: 'c-x', role: 'columnheader' }, el('span', { class: 'sr' }, 'Details')));
  const ths = cols.map(([k, l, cls]) => {
    const b = el('button', { type: 'button' }, l);
    const th = el('th', { class: cls, scope: 'col', role: 'columnheader', 'aria-sort': 'none' }, b);
    b.addEventListener('click', () => { st.dir = st.key === k ? -st.dir : (k === 'tic' || k === 'period' ? 1 : -1); st.key = k; sortSel.value = st.key + ':' + st.dir; draw(); });
    htr.append(th);
    return [k, th];
  });
  thead.append(htr); tab.append(thead, tbody);
  const count = document.getElementById('ccount');
  const more = el('button', { type: 'button', class: 'btn', style: 'margin-top:12px' }, 'Show more');
  tab.after(more);
  let limit = 60;
  more.addEventListener('click', () => { limit += 60; draw(); });

  // sort dropdown (shown in card mode, where there are no column headers)
  const SORTS = [['rank:-1', 'Most interesting first'], ['sde:-1', 'SDE, high → low'], ['period:1', 'Period, short → long'], ['period:-1', 'Period, long → short'],
    ['depth_ppm:-1', 'Depth, deep → shallow'], ['depth_ppm:1', 'Depth, shallow → deep'], ['rp_earth:1', 'Radius, small → large'], ['tmag:1', 'Brightest star first'], ['tic:1', 'TIC number']];
  const sortSel = el('select', { 'aria-label': 'Sort candidates by' }, SORTS.map(([v, l]) => el('option', { value: v }, l)));
  sortSel.value = `${st.key}:${st.dir}`;
  if (!sortSel.value) sortSel.value = 'rank:-1';
  sortSel.addEventListener('change', () => { const [k, dd] = sortSel.value.split(':'); st.key = k; st.dir = +dd; draw(); });
  F.prepend(el('div', { class: 'ctl sortctl' }, el('span', { class: 'lbl' }, 'Sort'), sortSel));

  const open = new Set();
  const detail = (c) => {
    const kv = [['Sector', c.sector], ['Match', c.known || (c.toi_host ? `TOI ${c.toi_host} host, other period` : 'none')], ['Duration', `${fmt(c.duration_h, 1)} h`],
      ['SNR', fmt(c.snr, 1)], ['Transits', c.n_transits], ['Tmag', fmt(c.tmag, 1)], ['Odd/even', `${fmt(c.odd_even_sigma, 1)}σ`], ['Phase 0.5', `${fmt(c.secondary_sigma, 1)}σ`],
      ['Exact period', `${fmt(c.period, 5)} d`], ['Exact depth', `${fmt(c.depth_ppm)} ppm`]];
    if (c.deep_vetted) kv.push(['Centroid offset', `${fmt(c.centroid_offset_px, 2)} px`], ['V-ratio', fmt(c.v_ratio, 2)], ['Size', `${fmt(c.rp_rjup, 2)} R_J`]);
    return el('div', { class: 'det' }, el('dl', {}, kv.map(([k, v]) => el('div', {}, el('dt', {}, k), el('dd', {}, String(v ?? '–'))))),
      el('a', { class: 'btn primary', href: `candidate.html#s${c.sector}/${c.tic}` }, 'Open light curve & checks →'),
      el('a', { class: 'btn', href: c.exofop, target: '_blank', rel: 'noopener' }, 'ExoFOP ↗'));
  };

  function draw() {
    const rows = all.filter((c) => (st.status === 'all' || c.status === st.status) && (!st.q || c.search.includes(st.q)) &&
      (!st.sector || c.sector === +st.sector) && (!st.flag || (st.flag === '-' ? !c.allflags : c.allflags.split(',').includes(st.flag))) &&
      (c.period == null || (c.period >= st.pmin && c.period <= st.pmax)) && (st.rmax >= 40 || c.rp_earth == null || c.rp_earth <= st.rmax) &&
      (st.tmax >= 17 || c.tmag == null || c.tmag <= st.tmax) && (c.sde ?? 0) >= st.sdemin)
      .sort((a, b) => { const x = a[st.key], y = b[st.key]; return (x == null) - (y == null) || (x > y ? 1 : x < y ? -1 : 0) * st.dir; });
    ths.forEach(([k, th]) => th.setAttribute('aria-sort', k === st.key ? (st.dir > 0 ? 'ascending' : 'descending') : 'none'));
    const out = [];
    for (const c of rows.slice(0, limit)) {
      const id = `${c.sector}-${c.tic}`, isOpen = open.has(id), did = `det-${id}`;
      const btn = el('button', { type: 'button', class: 'xbtn', 'aria-expanded': String(isOpen), 'aria-controls': did, 'aria-label': `More about TIC ${c.tic}` }, isOpen ? '−' : '+');
      btn.addEventListener('click', () => { isOpen ? open.delete(id) : open.add(id); draw(); });
      const cells = {
        tic: `<a href="candidate.html#s${c.sector}/${c.tic}">${c.tic}</a><small>S${c.sector}${c.known ? ' · ' + c.known : ''}</small>`,
        status: `<span class="badge b-${c.status}">${LABEL[c.status] || c.status}</span>`,
        period: per(c.period), depth_ppm: pct(c.depth_ppm), rp_earth: c.rp_earth == null ? '–' : c.rp_earth >= 10 ? fmt(c.rp_earth, 0) : fmt(c.rp_earth, 1),
        sde: fmt(c.sde, 1), allflags: flagChips(c),
      };
      out.push(el('tr', { role: 'row', class: isOpen ? 'open' : '' }, el('td', { role: 'cell', class: 'c-x' }, btn),
        cols.map(([k, l, cls]) => el('td', { role: 'cell', class: cls, 'data-label': l, html: `<span class="v">${cells[k]}</span>` }))));
      if (isOpen) out.push(el('tr', { role: 'row', class: 'detrow', id: did }, el('td', { role: 'cell', colspan: cols.length + 1 }, detail(c))));
    }
    tbody.replaceChildren(...out);
    more.hidden = rows.length <= limit;
    more.textContent = `Show more (${rows.length - Math.min(limit, rows.length)} left)`;
    count.textContent = `showing ${Math.min(limit, rows.length)} of ${rows.length} matching · ${all.length} listed`;
    save();
  }
  draw();
  return { setSector(s) { st.sector = String(s); secSel.value = st.sector; draw(); } };
}

function charts(d) {
  if (!window.Chart) return;
  const t = d.totals;
  new Chart(cRecovery, { type: 'doughnut', data: { labels: ['re-found & passed', 'found, failed a screen', 'missed'],
    datasets: [{ data: [t.toi_recovered, t.toi_period_found - t.toi_recovered, t.toi_hosts - t.toi_period_found], backgroundColor: [C.ok, C.warn, 'rgba(139,147,179,.35)'], borderWidth: 0 }] },
  options: { maintainAspectRatio: false, cutout: '68%', plugins: { legend: { position: 'right' } } } });
  const fu = d.deep_funnel;
  if (fu) {
    new Chart(cFunnel, { type: 'bar', data: { labels: ['passed screens, no match', 'deep-vetted', 'on-target (no blend)', 'not binary-like', 'survive everything'],
      datasets: [{ data: [fu.unmatched_passed, fu.deep_vetted, fu.survive_blend, fu.survive_eb, fu.survive_all],
        backgroundColor: ['rgba(124,196,255,.45)', 'rgba(124,196,255,.65)', C.acc, C.sub, C.acc2], borderWidth: 0 }] },
    options: { indexAxis: 'y', maintainAspectRatio: false, plugins: { legend: { display: false } } } });
    const de = Object.entries(fu.by_flag).filter((x) => x[1] > 0).sort((a, b) => b[1] - a[1]);
    new Chart(cDeep, { type: 'bar', data: { labels: de.map((x) => FLAG_TXT[x[0]] || x[0]), datasets: [{ data: de.map((x) => x[1]), backgroundColor: 'rgba(192,132,252,.6)', borderWidth: 0 }] },
      options: { indexAxis: 'y', maintainAspectRatio: false, plugins: { legend: { display: false } } } });
  }
  const sh = d.charts.sde_hist, ph = d.charts.period_hist;
  new Chart(cSde, { type: 'bar', data: { labels: sh.edges.slice(0, -1).map((x) => x.toFixed(0)), datasets: [{ data: sh.counts,
    backgroundColor: sh.edges.slice(0, -1).map((x) => (x >= 9 ? C.acc2 : 'rgba(124,196,255,.55)')), borderWidth: 0, barPercentage: 1, categoryPercentage: 1 }] },
  options: { maintainAspectRatio: false, plugins: { legend: { display: false } }, scales: { y: { type: 'logarithmic', title: { display: true, text: 'stars (log)' } }, x: { title: { display: true, text: 'SDE (orange = above threshold 9)' } } } } });
  new Chart(cPeriod, { type: 'bar', data: { labels: ph.edges.slice(0, -1).map((x) => x.toFixed(2)), datasets: [{ data: ph.counts, backgroundColor: 'rgba(192,132,252,.55)', borderWidth: 0, barPercentage: 1, categoryPercentage: 1 }] },
    options: { maintainAspectRatio: false, plugins: { legend: { display: false } }, scales: { x: { title: { display: true, text: 'period (days, log bins) — pile-ups are spacecraft systematics' } } } } });
  const fl = Object.entries(d.charts.flags).sort((a, b) => b[1] - a[1]);
  const FL = { low_sde: 'weak periodogram peak', too_deep: 'too deep (> 5%)', low_snr: 'low transit SNR', too_long: 'too long for the star', few_transits: '< 3 transits', deep: 'too deep (> 5%)',
    odd_even: 'odd ≠ even depth', variable: 'variable star', secondary: 'secondary eclipse', rp_too_large: 'too big for a planet', common_period: 'shared systematic period', no_dip: 'no dip' };
  new Chart(cFlags, { type: 'bar', data: { labels: fl.map((x) => FL[x[0]] || x[0]), datasets: [{ data: fl.map((x) => x[1]), backgroundColor: 'rgba(248,113,113,.55)', borderWidth: 0 }] },
    options: { indexAxis: 'y', maintainAspectRatio: false, plugins: { legend: { display: false } } } });
  new Chart(cSky, { type: 'scatter', data: { datasets: [{ data: d.charts.sky.map(([x, y]) => ({ x, y })), pointRadius: 1.2, backgroundColor: 'rgba(124,196,255,.55)' }] },
    options: { maintainAspectRatio: false, plugins: { legend: { display: false } }, scales: { x: { min: 0, max: 360, reverse: true, title: { display: true, text: 'RA (deg)' } }, y: { min: -90, max: 90, title: { display: true, text: 'Dec (deg)' } } } } });
  new Chart(cRuns, { type: 'bar', data: { labels: d.runs.map((r) => `S${r.sector} · ${r.started_at.slice(0, 10)}`), datasets: [
    { label: 'stars searched', data: d.runs.map((r) => r.n_searched), backgroundColor: 'rgba(124,196,255,.6)', yAxisID: 'y' },
    { type: 'line', label: 'median seconds per star', data: d.runs.map((r) => r.median_star_s), borderColor: C.acc2, backgroundColor: C.acc2, yAxisID: 'y1' }] },
  options: { maintainAspectRatio: false, scales: { y: { position: 'left' }, y1: { position: 'right', grid: { display: false } } } } });
}

function runs(d) {
  rtable(document.getElementById('runs'),
    [['a', 'Started (UTC)'], ['s', 'Sector'], ['n', 'Stars', 'num'], ['e', 'Errors', 'num'], ['w', 'Search', 'num'], ['m', 's / star', 'num'], ['b', 'Backend']],
    [...d.runs].reverse().map((r) => [r.started_at.slice(5, 16).replace('T', ' '), `S${r.sector}`, short(r.n_searched), r.n_errors,
      `${fmt(r.search_wall_s / 60, 0)} min`, fmt(r.median_star_s, 2), `${r.backend}${r.workers ? ` ×${r.workers}` : ''}`]));
}
