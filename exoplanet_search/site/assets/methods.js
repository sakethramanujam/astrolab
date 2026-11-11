// Methods & reliability: reads data/reliability.json (provisional binned summaries or the analysis PR's output).
import { el, fmt, slider, segmented, readout, surface, axes, palette, alpha, starfield } from './ui.js';

starfield();
const REPO = 'https://github.com/sakethramanujam/astrolab';
const GROUPS = { s70: 'Sector 70', s106_107: 'Sectors 106–107' };
const SPLITS = { tuning: 'Tuning stars', holdout: 'Held-out stars', all: 'All' };
const CHECK_INFO = {
  few_transits: 'fewer than 3 transits', odd_even: 'odd and even depths differ', secondary: 'dip at half phase',
  too_deep: 'deeper than 5%', too_long: 'too long for the star', variable: 'brightening at half phase',
  rp_too_large: 'bigger than 25 Earth radii', no_dip: 'no dip at all', common_period: 'period shared by many stars',
  centroid_offset: 'light moves off the star in transit', dip_off_target: 'dip sits on another pixel',
  known_eb: 'catalogued eclipsing binary', v_shape: 'V-shaped (grazing binary)', secondary_scan: 'faint secondary eclipse',
  rp_gt_2rjup: 'bigger than 2 Jupiters', density_mismatch: 'transit shape implies the wrong star',
  unresolved_neighbour: 'bright neighbour too close to rule out', centroid_shift_minor: 'small centroid shift',
  ruwe_high: 'Gaia RUWE > 1.4 (possible binary)', gaia_nss: 'Gaia non-single-star solution',
};
const BUCKETS = { few_transits: 'Fewer than 2 transits in the sector', long_duration: 'Transit too long for the 0.75-day detrending window',
  masked: 'Transits fell in flagged data or gap edges', below_threshold: 'Signal below our detection threshold',
  outcompeted: 'A stronger signal on the same star won', other: 'Other' };

const pct = (k, n) => (n ? `${Math.round((100 * k) / n)}%` : '–');
const J = await fetch('data/reliability.json', { cache: 'no-cache' }).then((r) => (r.ok ? r.json() : null)).catch(() => null);

// ---------------------------------------------------------------- header
const ver = document.getElementById('ver');
if (J) {
  ver.append(el('span', { class: 'pill' }, `pipeline v${J.pipeline_version}`), el('span', { class: 'pill' }, `config ${J.config_hash}`),
    el('a', { href: `${REPO}/tree/exosearch-v${J.pipeline_version}/exoplanet_search` }, 'frozen code'),
    el('a', { href: `${REPO}/blob/exoplanet-search/exoplanet_search/reproduce.md` }, 'reproduce'),
    el('a', { href: `${REPO}/blob/exoplanet-search/exoplanet_search/analysis/SCHEMA.md` }, 'data schema'));
  const DT = J.detection_threshold || {};
  if (DT.interim || DT.provisional) {
    const p = document.getElementById('prov');
    p.hidden = false;
    const secs = Object.entries(DT.by_sector || {});
    const ci = (a) => (a ? `${fmt(a[0], 1)}–${fmt(a[1], 1)}` : '–');
    const pop = DT.population || {};
    p.append(el('p', { html: DT.interim
      ? `<b>Interim, provisional threshold</b> (until both the flipped and the scrambled noise-only runs are complete). ` + secs.map(([k, v]) => `Sector ${k}: <b>S* = ${fmt(v.s_star, 1)}</b> (95% statistical interval ${ci(v.ci95)}` +
          (v.sys_halfwidth != null ? `; ±${fmt(v.sys_halfwidth, 1)} systematic from the tail starting point u, range ${ci(v.sys_range)} over u = ${(v.u_stable || []).join(', ')}` : '') + `), ` +
          `from a generalized-Pareto fit to the tail of the flipped light curves that pass every vetting check, ` +
          `${fmt(v.n_null)} tuning stars (${v.k_tail} vetted peaks above u = ${fmt(v.u, 2)}), solved for one false alarm per ${fmt(v.n_stars)} searched stars.`).join(' ') +
        ` False alarms are counted <i>after</i> vetting, and completeness below uses the same S* and the same vetting. ` +
        `A frozen star-quality cut, <b>registered after unblinding</b> (${DT.quality_cut}; see section 7 for its sensitivity), is applied to every population` +
        (Object.keys(pop).length ? ` (${Object.entries(pop).filter(([, v]) => v.measured).map(([k, v]) => `S${k}: ${fmt(v.dropped)} of ${fmt(v.searched)} dropped`).join(', ')})` : '') +
        `. The null run is still growing; the analysis code replaces these fits.`
      : `<b>Provisional.</b> The detection threshold is not calibrated yet.` }));
  }
} else {
  document.getElementById('main').prepend(el('div', { class: 'card' }, 'Reliability data not published yet.'));
}

function canvasFor(id, aspect, minH, maxH, draw) {
  const cv = document.getElementById(id);
  const s = surface(cv, { aspect, minH, maxH }, () => draw(s));
  return s;
}

// ---------------------------------------------------------------- completeness
const C = J?.completeness;
if (C) {
  const st = { outcome: 'passed', group: 's70', split: 'all', tb: 'all' };
  const groups = [...new Set(Object.keys(C.maps).map((k) => k.split('|')[0]))];
  st.group = groups.includes('s70') ? 's70' : groups[0];
  const ctl = document.getElementById('cctl');
  const tm = C.tmag_edges;
  const tbOpts = [{ value: 'all', label: 'All' }, ...tm.slice(0, -1).map((a, i) => ({ value: String(i),
    label: i === 0 ? `<${tm[1]}` : i === tm.length - 2 ? `>${tm[i]}` : `${tm[i]}–${tm[i + 1]}` }))];
  ctl.append(
    segmented({ id: 'co', label: 'Outcome', value: st.outcome, options: [{ value: 'detected', label: 'Detected' }, { value: 'passed', label: '+ passed vetting' }], onChange: (v) => { st.outcome = v; redraw(); } }).el,
    segmented({ id: 'cg', label: 'Sectors', value: st.group, options: groups.map((g) => ({ value: g, label: GROUPS[g] || g })), onChange: (v) => { st.group = v; redraw(); } }).el,
    segmented({ id: 'cs', label: 'Stars', value: st.split, options: Object.entries(SPLITS).map(([value, label]) => ({ value, label })), onChange: (v) => { st.split = v; redraw(); } }).el,
    segmented({ id: 'ct', label: 'TESS magnitude', value: st.tb, options: tbOpts, onChange: (v) => { st.tb = v; redraw(); } }).el);
  const ro = readout([{ k: 'cell', label: 'Cell' }, { k: 'n', label: 'Injected' }, { k: 'f', label: 'Recovered' }, { k: 'tot', label: 'Whole map' }]);
  document.getElementById('cread').append(ro.el);
  let hover = null;
  const P = palette();
  const map = () => C.maps[`${st.group}|${st.split}|${st.tb}`];
  const cell = { box: null };
  const s1 = canvasFor('cmap', 0.62, 260, 520, drawMap);
  const s2 = canvasFor('csnr', 0.42, 200, 340, drawSnr);
  function color(f) { // dark blue -> teal -> amber
    const a = [18, 30, 70], b = [94, 234, 212], c = [255, 184, 107];
    const m = f < 0.5 ? [a, b, f / 0.5] : [b, c, (f - 0.5) / 0.5];
    return `rgb(${m[0].map((v, i) => Math.round(v + (m[1][i] - v) * m[2])).join(',')})`;
  }
  function drawMap(s) {
    const { ctx, w, h } = s; ctx.clearRect(0, 0, w, h);
    const m = map(); if (!m) return;
    const box = { x0: 56, y0: 10, w: w - 70, h: h - 56 }; cell.box = box;
    const cw = box.w / 8, ch = box.h / 8;
    ctx.font = '11px "JetBrains Mono", monospace'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    for (let r = 0; r < 8; r++) for (let p = 0; p < 8; p++) {
      const n = m.n[r][p], k = m[st.outcome][r][p];
      const x = box.x0 + p * cw, y = box.y0 + (7 - r) * ch;
      ctx.fillStyle = n ? color(k / n) : 'rgba(255,255,255,.03)';
      ctx.fillRect(x + 1, y + 1, cw - 2, ch - 2);
      if (hover && hover[0] === r && hover[1] === p) { ctx.strokeStyle = P.fg; ctx.lineWidth = 2; ctx.strokeRect(x + 1, y + 1, cw - 2, ch - 2); }
      if (n && cw > 34) { ctx.fillStyle = k / n > 0.55 ? '#05070f' : P.fg; ctx.fillText(pct(k, n), x + cw / 2, y + ch / 2); }
    }
    ctx.fillStyle = P.muted; ctx.textBaseline = 'top';
    C.p_edges.forEach((v, i) => { if (i % 2 === 0 || cw > 44) ctx.fillText(v < 1 ? v.toFixed(1) : v.toFixed(0), box.x0 + i * cw, box.y0 + box.h + 4); });
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    C.r_edges.forEach((v, i) => ctx.fillText(v < 1 ? v.toFixed(1) : v.toFixed(0), box.x0 - 6, box.y0 + (8 - i) * ch));
    ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
    ctx.fillText('orbital period (days)', box.x0 + box.w / 2, h - 2);
    ctx.save(); ctx.translate(12, box.y0 + box.h / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'top'; ctx.fillText('planet radius (Earth = 1)', 0, 0); ctx.restore();
    let N = 0, K = 0; m.n.flat().forEach((v, i) => { N += v; K += m[st.outcome].flat()[i]; });
    ro.set('tot', `${pct(K, N)} of ${fmt(N)}`);
    if (hover) {
      const [r, p] = hover, n = m.n[r][p], k = m[st.outcome][r][p];
      ro.set('cell', `${C.p_edges[p].toFixed(1)}–${C.p_edges[p + 1].toFixed(1)} d · ${C.r_edges[r].toFixed(1)}–${C.r_edges[r + 1].toFixed(1)} R⊕`);
      ro.set('n', fmt(n)); ro.set('f', n ? `${pct(k, n)} (${k})` : '–');
    } else { ro.set('cell', 'hover or tap a cell'); ro.set('n', '–'); ro.set('f', '–'); }
  }
  function drawSnr(s) {
    const { ctx, w, h } = s; ctx.clearRect(0, 0, w, h);
    const cur = C.snr_curves[`${st.group}|${st.split}`]; if (!cur) return;
    const E = C.snr_edges, xmax = 40;
    const A = axes(ctx, { x0: 50, y0: 10, w: w - 64, h: h - 50 }, { x: [0, xmax], y: [0, 1], fy: (v) => `${Math.round(v * 100)}%`, xlabel: 'expected S of the injected planet', ylabel: 'recovered' });
    for (const [key, col] of [['det', P.acc], ['pas', P.acc2]]) {
      ctx.fillStyle = col; ctx.strokeStyle = alpha(col, 0.6);
      for (const [b, [n, det, pas]] of Object.entries(cur)) {
        const xc = (E[+b] + E[+b + 1]) / 2; if (xc > xmax || n < 3) continue;
        const k = key === 'det' ? det : pas, f = k / n, e = Math.sqrt(Math.max(f * (1 - f), 0.25 / n) / n);
        ctx.beginPath(); ctx.moveTo(A.X(xc), A.Y(Math.max(0, f - e))); ctx.lineTo(A.X(xc), A.Y(Math.min(1, f + e))); ctx.stroke();
        ctx.beginPath(); ctx.arc(A.X(xc), A.Y(f), 3.5, 0, 7); ctx.fill();
      }
    }
    const thr = J.detection_threshold?.default;
    if (thr) { ctx.strokeStyle = alpha(P.bad, 0.7); ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(A.X(thr), 10); ctx.lineTo(A.X(thr), h - 40); ctx.stroke(); ctx.setLineDash([]); }
    ctx.font = '12px Inter, sans-serif'; ctx.textAlign = 'right'; ctx.fillStyle = P.acc; ctx.fillText('● detected', w - 16, 22);
    ctx.fillStyle = P.acc2; ctx.fillText('● + passed vetting', w - 16, 40);
  }
  const pick = (e) => {
    const b = cell.box; if (!b) return;
    const r = s1.canvas.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
    const p = Math.floor((x - b.x0) / (b.w / 8)), rr = 7 - Math.floor((y - b.y0) / (b.h / 8));
    hover = p >= 0 && p < 8 && rr >= 0 && rr < 8 ? [rr, p] : null; drawMap(s1);
  };
  s1.canvas.addEventListener('pointermove', pick); s1.canvas.addEventListener('pointerdown', pick);
  function redraw() { drawMap(s1); drawSnr(s2); note(); }
  function note() {
    const sm = C.summary[`${st.group}|${st.split}`]; const n = document.getElementById('cnote'); n.replaceChildren();
    if (!sm) return;
    const fails = Object.entries(sm.fail_counts || {}).filter(([, v]) => v).sort((a, b) => b[1] - a[1]).slice(0, 5);
    n.append(el('p', { class: 'muted', html: `${fmt(sm.n)} injections in this selection: ${pct(sm.detected, sm.n)} detected, ${pct(sm.passed, sm.n)} also passed every check. ` +
      (Object.keys(sm.aliases || {}).length ? `Found at an alias period instead (not counted): ${Object.entries(sm.aliases).map(([k, v]) => `${k} ${v}`).join(', ')}. ` : '') +
      (fails.length ? `Checks that most often rejected detected fake planets: ${fails.map(([k, v]) => `<b>${k}</b> ${v}`).join(', ')}.` : '') }));
  }
  note();
} else document.getElementById('completeness').append(el('p', { class: 'muted' }, 'Injection runs in progress.'));

// ---------------------------------------------------------------- null tests
const N = J?.nulls;
if (N) {
  const keys = Object.keys(N.hist);
  const sectors = [...new Set(keys.map((k) => k.split('|')[0]))];
  const modes = [...new Set(keys.map((k) => k.split('|')[1]))];
  const st = { sector: sectors[0], mode: modes.includes('invert') ? 'invert' : modes[0], thr: J.detection_threshold?.by_sector?.[sectors[0]]?.s_star ?? J.detection_threshold?.default ?? 12 };
  const ctl = document.getElementById('nctl');
  const ro = readout([{ k: 'n', label: 'Stars' }, { k: 'above', label: 'Noise peaks ≥ S*' }, { k: 'vet', label: '…passing all checks' }, { k: 'fa', label: 'False alarms / sector' }]);
  document.getElementById('nread').append(ro.el);
  const P = palette();
  const s = canvasFor('nhist', 0.5, 220, 380, draw);
  ctl.append(
    segmented({ id: 'ns', label: 'Sector', value: st.sector, options: sectors.map((x) => ({ value: x, label: `S${x}` })), onChange: (v) => { st.sector = v; draw(s); } }).el,
    segmented({ id: 'nm', label: 'Noise model', value: st.mode, options: modes.map((m) => ({ value: m, label: m === 'invert' ? 'Flipped' : 'Scrambled' })), onChange: (v) => { st.mode = v; draw(s); } }).el,
    slider({ id: 'nt', label: 'Threshold S*', min: 5, max: 25, step: 0.25, value: st.thr, fmt: (v) => v.toFixed(2), onInput: (v) => { st.thr = v; draw(s); } }).el);
  function draw(sf) {
    const { ctx, w, h } = sf; ctx.clearRect(0, 0, w, h);
    const H = N.hist[`${st.sector}|${st.mode}`]; if (!H) return;
    const E = N.s_edges, xmax = 25, nr = H.n_real || 1;
    const ymax = Math.log10(Math.max(...H.all) + 1) + 0.2;
    const A = axes(ctx, { x0: 50, y0: 10, w: w - 64, h: h - 50 }, { x: [0, xmax], y: [0, ymax], fy: (v) => (Number.isInteger(v) ? `10${'⁰¹²³⁴⁵'[v] ?? v}` : ''), xlabel: 'best-peak S of each star', ylabel: 'stars' });
    const bw = A.X(E[1]) - A.X(E[0]);
    H.all.forEach((c, i) => { if (!c || E[i] >= xmax) return; ctx.fillStyle = E[i] >= st.thr ? alpha(P.bad, 0.75) : alpha(P.acc, 0.45); const y = A.Y(Math.log10(c + 1)); ctx.fillRect(A.X(E[i]), y, Math.max(1, bw - 0.5), A.Y(0) - y); });
    H.vetted.forEach((c, i) => { if (!c || E[i] >= xmax) return; ctx.fillStyle = P.acc2; const y = A.Y(Math.log10(c + 1)); ctx.fillRect(A.X(E[i]) + bw * 0.25, y, Math.max(1, bw * 0.5), A.Y(0) - y); });
    ctx.strokeStyle = P.bad; ctx.beginPath(); ctx.moveTo(A.X(st.thr), 10); ctx.lineTo(A.X(st.thr), h - 40); ctx.stroke();
    let above = 0, vet = 0; H.all.forEach((c, i) => { if (E[i] >= st.thr) above += c; }); H.vetted.forEach((c, i) => { if (E[i] >= st.thr) vet += c; });
    const n = N.n_stars[`${st.sector}|${st.mode}`];
    const nsec = N.threshold?.[st.sector]?.n_stars;
    ro.set('n', fmt(n)); ro.set('above', fmt(above / nr, nr > 1 ? 1 : 0)); ro.set('vet', fmt(vet / nr, nr > 1 ? 1 : 0), vet ? 'warnv' : 'good');
    const fa = nsec && n ? (vet / nr) * (nsec / n) : vet / nr;      // scaled from the null sample to the whole sector
    ro.set('fa', vet ? fmt(fa, 1) : `< ${fmt(nsec && n ? nsec / n : 1, 1)}`, fa > 1 ? 'badv' : vet ? 'warnv' : 'good');
    ctx.font = '12px Inter, sans-serif'; ctx.textAlign = 'right'; ctx.fillStyle = P.acc2; ctx.fillText('▮ passes every other check', w - 16, 22);
  }
  const t = document.getElementById('nthr');
  const rows = Object.entries(N.threshold || {}).filter(([k]) => !k.includes(':') && N.threshold[k].ci68).sort();
  const ci = (a) => (a ? `${fmt(a[0], 1)}–${fmt(a[1], 1)}` : '–');
  for (const [k, v] of rows) {
    const sb = N.threshold_by_null?.[k];
    const scr = sb?.scramble?.s_star ? `${fmt(sb.scramble.s_star, 1)} (95% ${ci(sb.scramble.ci95)}, ${fmt(sb.scramble.n)} stars)` : 'pending (run not finished)';
    t.append(el('p', { class: 'muted', html: `<b>Provisional interim S* for sector ${k}: ${fmt(v.s_star, 1)}</b> (95% statistical ${ci(v.ci95)}${v.sys_halfwidth != null ? `, ±${fmt(v.sys_halfwidth, 1)} systematic` : ''}). ` +
      `Side by side: flipped null ${fmt(v.s_star, 1)}; scrambled null ${scr}. It stays provisional until both are complete. ` +
      `Fit: generalized Pareto above S = ${fmt(v.u, 2)} to ${v.k_tail} vetted noise peaks among ${fmt(v.n)} tuning stars (shape ξ = ${fmt(v.xi, 2)}); ` +
      `the largest vetted noise peak so far is S = ${fmt(v.max_vetted, 1)}. Cameras are pooled. The histogram shows stars that pass the quality cut; the slider starts at S*.` }));
  }
} else document.getElementById('falsealarms').append(el('p', { class: 'muted' }, 'Null runs in progress.'));

// ---------------------------------------------------------------- per-check rates
const K = J?.check_rates;
if (K) {
  const keys = Object.keys(K);
  const st = { group: keys.some((k) => k.startsWith('s70')) ? 's70' : keys[0].split('|')[0], split: 'tuning' };
  const ctl = document.getElementById('kctl');
  ctl.append(
    segmented({ id: 'kg', label: 'Sectors', value: st.group, options: [...new Set(keys.map((k) => k.split('|')[0]))].map((g) => ({ value: g, label: GROUPS[g] || g })), onChange: (v) => { st.group = v; draw(); } }).el,
    segmented({ id: 'ks', label: 'Stars', value: st.split, options: [{ value: 'tuning', label: 'Tuning' }, { value: 'holdout', label: 'Held-out' }], onChange: (v) => { st.split = v; draw(); } }).el);
  const list = document.getElementById('klist');
  const bar = (r, cls) => { const f = r.n ? r.k / r.n : 0; const [lo, hi] = r.ci68 || [null, null];
    return el('div', { class: `bar ${cls}` }, el('i', { style: `width:${(f * 100).toFixed(1)}%` }),
      lo != null ? el('em', { style: `left:${(lo * 100).toFixed(1)}%;width:${((hi - lo) * 100).toFixed(1)}%` }) : null,
      el('span', {}, `${pct(r.k, r.n)} · ${r.k}/${r.n}`)); };
  function draw() {
    list.replaceChildren();
    const R = K[`${st.group}|${st.split}`];
    if (!R) { list.append(el('p', { class: 'muted' }, 'No labelled signals in this selection yet.')); return; }
    list.append(el('div', { class: 'brow head' }, el('b', {}, 'Check'), el('b', {}, 'Fires on planets (want ≈0)'), el('b', {}, 'Fires on false positives (want high)')));
    const order = Object.entries(R).filter(([k]) => k !== '_any_reject').sort((a, b) => (b[1].fp.n ? b[1].fp.k / b[1].fp.n : 0) - (a[1].fp.n ? a[1].fp.k / a[1].fp.n : 0));
    for (const [c, r] of [['_any_reject', { ...R._any_reject, reject: true }], ...order]) {
      list.append(el('div', { class: 'brow' },
        el('div', {}, el('b', {}, c === '_any_reject' ? 'Any rejecting check' : c), el('small', {}, c === '_any_reject' ? 'the whole chain' : `${CHECK_INFO[c] || ''}${r.reject ? '' : ' · warning only'}`)),
        bar(r.planet, 'pl'), bar(r.fp, 'fpb')));
    }
    list.append(el('p', { class: 'muted' }, 'Bars show the fraction; the lighter band is a 68% Wilson interval. Catalogue eclipsing binaries are left out of the known-EB row, since that check uses the same catalogue.'));
  }
  draw();
} else document.getElementById('checks').append(el('p', { class: 'muted' }, 'Labelled-signal checks in progress.'));

// ---------------------------------------------------------------- missed TOIs + survivors
const M = J?.missed_tois;
const ml = document.getElementById('mlist');
for (const [sec, m] of Object.entries(M || {})) {
  ml.append(el('p', {}, `Sector ${+sec.replace('s', '')}: ${m.n} known TOI hosts where our best peak was not the TOI's period (first matching reason, in this order):`));
  for (const [b, label] of Object.entries(BUCKETS)) {
    const k = m.buckets?.[b] || 0;
    ml.append(el('div', { class: 'brow two' }, el('div', {}, el('b', {}, label)), el('div', { class: 'bar fpb' }, el('i', { style: `width:${((100 * k) / m.n).toFixed(1)}%` }), el('span', {}, `${k} · ${pct(k, m.n)}`))));
  }
}
if (!M || !Object.keys(M).length) ml.append(el('p', { class: 'muted' }, 'Diagnosis in progress.'));

const S = J?.survivors;
const sl = document.getElementById('slist');
if (S) {
  sl.append(el('p', {}, S.note || ''));
  const g = el('div', { class: 'grid two' });
  for (const c of S.rows) {
    const tone = c.status === 'rejected' ? 'bad' : c.status.startsWith('marginal') ? 'warn' : '';
    g.append(el('a', { class: 'card scard', href: `candidate.html#s${c.sector}/${c.tic}` },
      el('div', { class: 'srow' }, el('b', {}, `TIC ${c.tic}`), el('span', { class: `pill ${tone}` }, c.status)),
      el('div', { class: 'muted' }, `P ${fmt(c.period, 3)} d · ${c.rp_earth ? `${fmt(c.rp_earth, 1)} R⊕` : 'size ?'} · S ${fmt(c.S, 1)}${c.s_star ? ` (S* ${fmt(c.s_star, 1)})` : ''} · FAP ${c.fap != null ? (c.fap < 1e-3 ? c.fap.toExponential(1) : fmt(c.fap, 3)) : '–'}`),
      c.expected_false_alarms_per_sector != null ? el('div', { class: 'muted' }, `Expected false alarms per sector at this S (after vetting): ${fmt(c.expected_false_alarms_per_sector, 2)}`) : null,
      el('div', { class: 'muted' }, `ρ★ transit/TIC ${c.rho_transit != null && c.rho_tic ? fmt(c.rho_transit / c.rho_tic, 2) : '–'} (z ${fmt(c.rho_z, 1)}) · RUWE ${fmt(c.ruwe, 2)}`),
      c.flags?.length ? el('div', { class: 'fchips' }, c.flags.map((f) => el('span', { class: 'fchip' }, f))) : null));
  }
  sl.append(g);
} else sl.append(el('p', { class: 'muted' }, 'Re-vetting in progress.'));

// ---------------------------------------------------------------- interim review (section 7)
const RV = J?.review, WN = J?.whitenoise_neff, FF = J?.ffi_check_318707987;
const rvb = document.getElementById('rvbody');
if (rvb && RV) {
  const tbl = (cap, head, rows) => el('table', { class: 'mini' }, el('caption', {}, cap),
    el('tr', {}, head.map((h) => el('th', {}, h))), rows.map((r) => el('tr', {}, r.map((x) => el('td', {}, x == null ? '–' : String(x))))));
  const g = RV.gpd || {};
  const ci = (a) => (a ? `${fmt(a[0], 1)}–${fmt(a[1], 1)}` : '–');
  rvb.append(el('h3', {}, 'Tail starting point'), el('p', { html:
    `The GPD starting point u is chosen by a fixed rule: ${g.u_rule}. Vetted peaks exist only at S ≥ 6, so u starts there. ` +
    `Chosen u = ${fmt(g.u, 2)} (${g.k} tail points): <b>S* = ${fmt(g.s_star, 1)}</b>, 95% statistical interval ${ci(g.ci95_stat)}, ` +
    `systematic ±${fmt(g.sys_halfwidth, 1)} (range ${ci(g.sys_range)} over u = ${(g.stable_us || []).join(', ')}). ` +
    `<a href="data/review/gpd_threshold_choice_s0070.png">Mean-residual-life and stability plots (PNG)</a>.` }));
  rvb.append(el('p', { html: `Cameras are pooled: each camera's tail is too thin (tens of points) for a stable fit, so no per-camera threshold is used or shown until the samples grow.` }));
  const qs = Object.entries(RV.quality_sensitivity || {});
  rvb.append(el('h3', {}, 'Star-quality cut: frozen, registered after unblinding'), el('p', { html:
    `The cut was written down after we had seen summary statistics of the flipped null, but before evaluating it on any population. It is frozen. Its effect at other thresholds:` }),
    tbl('Sensitivity to the β_max cut (sector 70, tuning stars, allowed radius range)', ['β_max cut', 'Dropped', 'S*', 'Vetted, S_exp ≥ 25'],
      qs.map(([k, v]) => [k === 'None' ? 'none' : k, fmt(v.stars_dropped), fmt(v.s_star, 1), `${fmt(100 * v.vetted_completeness_allowed_snr_ge_25, 0)}%`])));
  const ts = Object.entries(RV.toi_dropped_vs_kept || {});
  rvb.append(tbl('Known planets (TOI hosts) among dropped and kept stars', ['Sector', 'Dropped: hosts, re-found', 'Kept: hosts, re-found', 'Odds ratio (p)'],
    ts.map(([s, v]) => [s, `${v.dropped.toi_hosts} (${fmt(100 * v.dropped.toi_host_rate, 2)}%), ${v.dropped.recovered}`, `${v.kept.toi_hosts} (${fmt(100 * v.kept.toi_host_rate, 2)}%), ${v.kept.recovered}`,
      `${fmt(v.fisher_overrep.odds_ratio, 2)} (${v.fisher_overrep.p.toExponential(0)})`])),
    el('p', { class: 'muted' }, 'Known planets are under-represented among dropped stars, not over-represented.'));
  const ca = RV.completeness_allowed?.tuning || {};
  const hs = RV.high_snr_failures || {};
  rvb.append(el('h3', {}, 'Completeness inside the allowed size range'), el('p', { html:
    `The too_deep (depth ≥ 5%) and size (R ≥ 25 R⊕) caps reject large bodies by design, so completeness is quoted only for injections they allow ` +
    `(${fmt(100 * RV.allowed_fraction, 0)}% of injections). Of ${hs.n} high-S_exp (≥ 50) injections that were found but failed vetting, ` +
    `too_deep is involved in ${fmt(100 * hs.too_deep_any, 0)}% (the only failure in ${fmt(100 * hs.too_deep_only, 0)}%), odd_even in ${fmt(100 * hs.odd_even_any, 0)}%.` }),
    tbl('Sector 70, tuning stars, at the interim S*', ['Expected S', 'n', 'Found', 'Found and vetted'],
      Object.entries(ca).map(([k, v]) => [k.replace(/[\[\)]/g, '').replace(', ', '–').replace('1000000000.0', '∞'), v.n, `${fmt(100 * v.detected, 0)}%`, `${fmt(100 * v.vetted, 0)}%`])));
  const tf = RV.tail_fap || {};
  rvb.append(el('h3', {}, 'TIC 318707987: marginal'), el('p', { html:
    `At its S of 14.8 the tail fit gives ${fmt(tf['14.8']?.expected_false_alarms_per_sector, 2)} expected false alarms per sector after vetting ` +
    `(per-star ${tf['14.8']?.per_star?.toExponential(1)}), so it is <b>marginal</b>, not a detection. ` +
    (FF ? `The planned independent test in the sector 42/43 full-frame images is not possible: the star fell off the science pixels in both sectors ` +
      `(the image cutouts hold ~0.02 e/s where the star would give ~75,000 e/s; the TESS pointing tool lists sector 70 only). ` +
      `As a photometry cross-check only, our own aperture photometry of the sector 70 full-frame images finds the same period blind ` +
      `(P ${fmt(FF.control_blind?.period, 4)} d, S ${fmt(FF.control_blind?.S, 1)}, ${fmt(FF.control_blind?.depth_ppm, 0)} ppm), but those are the same photons. ` : '') +
    `<a href="${REPO}/blob/exoplanet-search/exoplanet_search/docs/preregistration/tic318707987_ffi.md">Test registration</a>.` }));
  if (WN) {
    const f = WN.fits?.S_oot || {};
    rvb.append(el('h3', {}, 'White noise: effective number of trials'), el('p', { html:
      `Pure white noise with the real sector-70 cadence, through the full chain, with the out-of-transit error term restored: Gumbel fit to the per-star maxima gives ` +
      `an expected maximum of ${fmt(f.expected_max, 2)} and about ${fmt(f.joint_N_eff / 1e4, 1)} × 10⁴ effective independent trials per star (per-trial spread ${fmt(f.joint_trial_sd, 2)}). ` +
      `The white-noise threshold for one false alarm per sector would be ${fmt(f.s_star_white_joint, 1)}; the real (vetted) S* is far higher, so the false-alarm tail comes from non-white residuals.` }));
  }
} else if (rvb) rvb.append(el('p', { class: 'muted' }, 'Review in progress.'));
