import { el, fmt, slider, toggle, readout, surface, axes, palette, alpha, notice, starfield } from './ui.js';

starfield();
const LABEL = { needs_review: 'Needs review', uncalibrated: 'Not calibrated', known: 'Known', submitted: 'Submitted', rejected: 'Rejected' };
const ppmFmt = (v) => (Math.abs(v) >= 10000 ? `${(v / 1e4).toFixed(Math.abs(v) >= 1e5 ? 0 : 1)}%` : `${Math.round(v)}`);

function parseHash() {
  const m = location.hash.match(/s(\d+)\/(\d+)/);
  return m ? { sector: +m[1], tic: +m[2] } : null;
}

async function main() {
  const d = await (await fetch('data/progress.json', { cache: 'no-store' })).json();
  const RANK = { submitted: 5, needs_review: 4, uncalibrated: 3, known: 2, rejected: 1 };
  const list = [...d.candidates].sort((a, b) => RANK[b.status] * 1000 + b.sde - (RANK[a.status] * 1000 + a.sde));
  const want = parseHash();
  const i = want ? list.findIndex((c) => c.tic === want.tic && c.sector === want.sector) : -1;
  if (i < 0) {
    document.getElementById('title').textContent = 'Candidate not found';
    document.getElementById('lede').innerHTML = 'It may have been re-classified or dropped from the public list. <a href="./#candidates">Back to all candidates</a>.';
    return;
  }
  render(list[i], list[i - 1], list[i + 1]);
}

function render(c, prev, next) {
  const name = `TIC ${c.tic}`;
  document.title = `${name} · S${c.sector} · Transit Hunt`;
  document.getElementById('crumb').textContent = `${name} · sector ${c.sector}`;
  document.getElementById('title').textContent = `${name}`;
  document.getElementById('lede').innerHTML = `<span class="badge b-${c.status}">${LABEL[c.status] || c.status}</span> ` +
    (c.known ? `Our blind search re-found <b>${c.known}</b>. ` : c.toi_host ? `This star hosts TOI ${c.toi_host}, but this signal is at a different period. ` : 'No catalogue match at this period. ') +
    `A ${fmt(c.depth_ppm)} ppm dip lasting ${fmt(c.duration_h, 1)} h, every ${fmt(c.period, 4)} days, seen ${c.n_transits} times in sector ${c.sector}.` +
    (c.status === 'needs_review' ? ` <b>Not a detection:</b> it passed our automated checks with S ${fmt(c.S, 1)} ≥ interim S* ${fmt(c.s_star, 1)}, which is where human vetting starts.` : '') +
    (c.status === 'uncalibrated' ? ` <b>Not calibrated:</b> ${c.calib_note || 'no interim S* comparison'}. It passed the v1.0.0 screens only.` : '');
  const mast = `https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html?searchQuery=${encodeURIComponent(`TIC ${c.tic}`)}`;
  document.getElementById('links').append(
    el('a', { class: 'btn', href: c.exofop, target: '_blank', rel: 'noopener' }, 'ExoFOP ↗'),
    el('a', { class: 'btn', href: mast, target: '_blank', rel: 'noopener' }, 'MAST portal ↗'),
    c.plot ? el('a', { class: 'btn', href: c.plot, target: '_blank', rel: 'noopener' }, 'Diagnostic sheet (PNG) ↗') : null,
    c.lc ? el('a', { class: 'btn', href: c.lc }, 'Light-curve data (JSON)') : null);

  const body = document.getElementById('body');
  const ro = readout([{ k: 'p', label: 'Period' }, { k: 'd', label: 'Depth' }, { k: 'dur', label: 'Duration' }, { k: 'rp', label: 'Radius' },
    { k: 'sde', label: 'SDE' }, { k: 'snr', label: 'SNR' }, { k: 'tmag', label: 'TESS mag' }, { k: 'star', label: 'Star' }]);
  ro.set('p', `${fmt(c.period, 5)} d`); ro.set('d', `${fmt(c.depth_ppm)} ppm`); ro.set('dur', `${fmt(c.duration_h, 2)} h`);
  ro.set('rp', c.rp_earth != null ? `${fmt(c.rp_earth, 1)} R⊕` : '–', c.rp_rjup > 2 ? 'badv' : '');
  ro.set('sde', fmt(c.sde, 1), c.sde >= 9 ? 'good' : 'badv'); ro.set('snr', fmt(c.snr, 1), c.snr >= 7 ? 'good' : 'badv');
  ro.set('tmag', fmt(c.tmag, 2)); ro.set('star', c.radius_rsun ? `${fmt(c.radius_rsun, 2)} R☉${c.teff ? ` · ${fmt(c.teff)} K` : ''}` : '–');
  body.append(el('div', { class: 'card' }, ro.el));

  const plots = el('div', { class: 'card' });
  body.append(plots);
  if (c.lc) fetch(c.lc).then((r) => r.json()).then((lc) => interactive(plots, c, lc)).catch(() => fallback(plots, c));
  else fallback(plots, c);

  body.append(el('div', { class: 'card', style: 'margin-top:16px' }, el('h2', {}, 'Checks'),
    el('p', { class: 'kbd', style: 'margin-top:-6px' }, 'Every automated screen and deep-vetting test this signal went through. Each one is a reason a dip could be something other than a planet.'),
    el('div', { class: 'checks' }, checks(c).map(check))));

  const np = document.getElementById('np');
  np.append(prev ? el('a', { href: `#s${prev.sector}/${prev.tic}` }, el('small', {}, '← Previous'), `TIC ${prev.tic} · ${LABEL[prev.status]}`) : el('a', { href: './#candidates' }, el('small', {}, '← Back'), 'All candidates'),
    next ? el('a', { class: 'n', href: `#s${next.sector}/${next.tic}` }, el('small', {}, 'Next →'), `TIC ${next.tic} · ${LABEL[next.status]}`) : el('a', { class: 'n', href: './#candidates' }, el('small', {}, 'Back →'), 'All candidates'));
}
addEventListener('hashchange', () => location.reload());

function fallback(root, c) {
  root.append(el('h2', {}, 'Diagnostic sheet'),
    c.plot ? el('img', { src: c.plot, alt: `Diagnostic plots for TIC ${c.tic}`, style: 'width:100%;border-radius:12px;background:#fff' }) : el('p', {}, 'No light-curve data published for this signal.'));
}

/* ───── interactive light-curve panels ───── */
function interactive(root, c, lc) {
  const st = { win: Math.min(Math.max(3 * lc.duration_h, 3), maxWin(lc)), bin: Math.max(10, Math.round(lc.duration_h * 60 / 6 / 5) * 5), oddeven: false, model: true, marks: true };
  const depth = lc.depth_ppm, dh = lc.duration_h;
  const allY = lc.near.ppm.concat(lc.sec.ppm).filter((v) => v != null).sort((a, b) => a - b);
  const q = (p) => allY[Math.floor(p * (allY.length - 1))] ?? 0;
  const yr = [Math.min(-1.5 * depth, q(0.005)), Math.max(0.6 * depth, q(0.995))];

  const mk = (title, legend, aspect = 0.42) => {
    const cv = el('canvas', { role: 'img', 'aria-label': title });
    const tip = el('div', { class: 'tip', 'aria-hidden': 'true' });
    root.append(el('div', { class: 'cap' }, el('b', {}, title), legend || ''), el('div', { class: 'plotwrap' }, cv, tip));
    return { cv, tip, s: surface(cv, { aspect, minH: 220, maxH: 380 }, () => draw()) };
  };
  const leg = (items) => el('span', { class: 'legend-inline' }, items.map(([c_, t]) => el('span', {}, el('i', { style: `background:${c_}` }), t)));
  const P = palette();

  const stage = el('div', { class: 'stage' });
  const viz = el('div', { class: 'viz' }), controls = el('aside', { class: 'controls', 'aria-label': 'Plot controls' });
  stage.append(viz, controls);
  root.append(el('h2', {}, 'Light curve, interactively'), stage);
  const saveRoot = root; root = viz;
  const F = mk('Folded on the period, around the transit', leg([[P.acc, 'each 10-min point'], [P.acc2, 'binned'], [P.ok, 'box model'], [P.sub, 'odd transits (if split)']]));
  const S = mk('Half an orbit later (phase 0.5): a secondary eclipse would mean a binary star', leg([[P.acc, 'points'], [P.acc2, 'binned'], [P.bad, 'depth if it were a twin star']]), 0.3);
  const T = mk('The whole sector, hourly bins', leg([[P.acc, 'brightness'], [P.acc2, 'predicted transits']]), 0.3);
  root = saveRoot;

  const sW = slider({ id: 'win', label: 'Fold window', min: Math.max(1, dh), max: maxWin(lc), step: 0.25, value: st.win, fmt: (v) => `±${v.toFixed(2)} h`, onInput: (v) => { st.win = v; draw(); },
    hint: 'Zoom in to judge the shape: a flat bottom suits a planet, a V suits two grazing stars.' });
  const sB = slider({ id: 'bin', label: 'Bin width', min: 10, max: 120, step: 5, value: st.bin, fmt: (v) => `${v} min`, onInput: (v) => { st.bin = v; draw(); },
    hint: 'Averaging neighbouring points beats down noise by √N.' });
  const tO = toggle({ id: 'oe', label: 'Split odd and even transits', value: false, onChange: (v) => { st.oddeven = v; draw(); } });
  const tM = toggle({ id: 'model', label: 'Show box model (BLS fit)', value: true, onChange: (v) => { st.model = v; draw(); } });
  const tK = toggle({ id: 'marks', label: 'Mark predicted transits', value: true, onChange: (v) => { st.marks = v; draw(); } });
  const roF = readout([{ k: 'in', label: 'In-transit mean' }, { k: 'odd', label: 'Odd depth' }, { k: 'even', label: 'Even depth' }, { k: 'sec', label: 'Phase-0.5 mean' }]);
  controls.append(el('h3', {}, 'Fold'), sW.el, sB.el, tO.el, tM.el, el('h3', {}, 'Sector view'), tK.el, el('h3', {}, 'Measured'), roF.el);
  root.append(notice([
    '<b>A planet gives a flat-bottomed U</b> once the whole disc is in front of the star; two grazing stars give a sharper V.',
    '<b>Odd and even transits should match.</b> If every other dip is deeper, the “planet” is probably two stars of different size eclipsing in turn, at twice this period.',
    '<b>Phase 0.5 should be flat.</b> A dip there is the second star being eclipsed; a planet is far too faint to show one in TESS data.',
    'In the sector view, the predicted transits should land on real dips, not on gaps or the ramps at the edges of the gaps.',
  ]));

  const box = (s) => ({ x0: 58, y0: 10, w: s.w - 70, h: s.h - 48 });
  const binned = (xs, ys, w, lo, hi) => {
    const n = Math.max(1, Math.round((hi - lo) / w));
    const sum = new Float64Array(n), cnt = new Uint32Array(n);
    xs.forEach((x, k) => { const y = ys[k]; if (y == null || x < lo || x >= hi) return; const j = Math.floor((x - lo) / w); sum[j] += y; cnt[j]++; });
    return [...sum].map((v, j) => [lo + (j + 0.5) * w, cnt[j] ? v / cnt[j] : null]);
  };
  const mean = (a) => (a.length ? a.reduce((x, y) => x + y, 0) / a.length : NaN);
  const hover = (P_, pts, X, Y, fmtPt) => {
    P_.cv.onpointermove = (e) => {
      const r = P_.cv.getBoundingClientRect(), mx = e.clientX - r.left, my = e.clientY - r.top;
      let best = null, bd = 400;
      for (const p of pts) { const dd = (X(p[0]) - mx) ** 2 + (Y(p[1]) - my) ** 2; if (dd < bd) { bd = dd; best = p; } }
      if (!best) { P_.tip.style.display = 'none'; return; }
      P_.tip.style.display = 'block'; P_.tip.style.left = X(best[0]) + 'px'; P_.tip.style.top = Y(best[1]) + 'px';
      P_.tip.textContent = fmtPt(best);
    };
    P_.cv.onpointerleave = () => { P_.tip.style.display = 'none'; };
  };

  function foldPanel(Pn, xs, ys, es, isSec) {
    const { ctx, w, h } = Pn.s;
    if (!w) return;
    ctx.clearRect(0, 0, w, h);
    const b = box(Pn.s);
    const { X, Y } = axes(ctx, b, { x: [-st.win, st.win], y: yr, xlabel: isSec ? 'hours from phase 0.5' : 'hours from mid-transit', ylabel: 'ppm', fy: ppmFmt });
    ctx.save(); ctx.beginPath(); ctx.rect(b.x0, b.y0, b.w, b.h); ctx.clip();
    const pts = [];
    xs.forEach((x, k) => {
      const y = ys[k];
      if (y == null || Math.abs(x) > st.win) return;
      const odd = es && es[k] % 2 !== 0;
      ctx.fillStyle = alpha(st.oddeven && es ? (odd ? P.sub : P.acc) : P.acc, 0.45);
      ctx.beginPath(); ctx.arc(X(x), Y(y), 1.8, 0, 6.283); ctx.fill();
      pts.push([x, y, es ? es[k] : null]);
    });
    const w_ = st.bin / 60;
    const line = (bb, col) => { ctx.strokeStyle = col; ctx.lineWidth = 2.4; ctx.beginPath(); let on = false;
      for (const [x, y] of bb) { if (y == null) { on = false; continue; } on ? ctx.lineTo(X(x), Y(y)) : ctx.moveTo(X(x), Y(y)); on = true; } ctx.stroke(); };
    if (st.oddeven && es) {
      for (const [par, col] of [[0, P.acc2], [1, P.sub]]) {
        const ii = xs.map((_, k) => k).filter((k) => Math.abs(es[k] % 2) === par);
        line(binned(ii.map((k) => xs[k]), ii.map((k) => ys[k]), w_, -st.win, st.win), col);
      }
    } else line(binned(xs, ys, w_, -st.win, st.win), P.acc2);
    if (isSec) {
      ctx.setLineDash([5, 4]); ctx.strokeStyle = P.bad; ctx.lineWidth = 1.2;
      ctx.beginPath(); ctx.moveTo(X(-dh / 2), Y(-depth)); ctx.lineTo(X(dh / 2), Y(-depth)); ctx.stroke(); ctx.setLineDash([]);
    } else if (st.model) {
      ctx.strokeStyle = P.ok; ctx.lineWidth = 1.6; ctx.beginPath();
      ctx.moveTo(X(-st.win), Y(0)); ctx.lineTo(X(-dh / 2), Y(0)); ctx.lineTo(X(-dh / 2), Y(-depth)); ctx.lineTo(X(dh / 2), Y(-depth)); ctx.lineTo(X(dh / 2), Y(0)); ctx.lineTo(X(st.win), Y(0)); ctx.stroke();
    }
    ctx.restore();
    hover(Pn, pts, X, Y, (p) => `${p[0].toFixed(2)} h · ${Math.round(p[1])} ppm${p[2] != null ? ` · transit #${p[2]}` : ''}`);
  }

  function seriesPanel() {
    const { ctx, w, h } = T.s;
    if (!w) return;
    ctx.clearRect(0, 0, w, h);
    const b = box(T.s);
    const ts = lc.series.t, ys = lc.series.ppm;
    const tmax = ts[ts.length - 1] || 27;
    const sy = ys.filter((v) => v != null).sort((a, b_) => a - b_);
    const lo = Math.min(sy[Math.floor(0.002 * sy.length)] ?? -depth, -depth * 1.3), hi = Math.max(sy[Math.floor(0.998 * sy.length)] ?? depth, depth * 0.4);
    const { X, Y } = axes(ctx, b, { x: [0, tmax], y: [lo, hi], xlabel: 'days since the first good point', ylabel: 'ppm', fy: ppmFmt });
    ctx.save(); ctx.beginPath(); ctx.rect(b.x0, b.y0, b.w, b.h); ctx.clip();
    if (st.marks) {
      ctx.strokeStyle = alpha(P.acc2, 0.5); ctx.lineWidth = 1;
      for (let t = lc.t0_offset - Math.ceil(lc.t0_offset / lc.period) * lc.period; t <= tmax; t += lc.period) {
        if (t < 0) continue; ctx.beginPath(); ctx.moveTo(X(t), b.y0); ctx.lineTo(X(t), b.y0 + 8); ctx.moveTo(X(t), b.y0 + b.h - 8); ctx.lineTo(X(t), b.y0 + b.h); ctx.stroke();
      }
    }
    ctx.strokeStyle = P.acc; ctx.lineWidth = 1.2; ctx.beginPath();
    let last = -1;
    const pts = [];
    ts.forEach((t, k) => { const y = ys[k]; if (y == null) return; (t - last > 0.1 || last < 0) ? ctx.moveTo(X(t), Y(y)) : ctx.lineTo(X(t), Y(y)); last = t; pts.push([t, y]); });
    ctx.stroke(); ctx.restore();
    hover(T, pts, X, Y, (p) => `day ${p[0].toFixed(2)} · ${Math.round(p[1])} ppm`);
  }

  function draw() {
    foldPanel(F, lc.near.h, lc.near.ppm, lc.near.e, false);
    foldPanel(S, lc.sec.h, lc.sec.ppm, null, true);
    seriesPanel();
    const inT = (xs, ys, f = () => true) => mean(xs.map((x, k) => (Math.abs(x) < dh / 2 && ys[k] != null && f(k) ? ys[k] : null)).filter((v) => v != null));
    const odd = inT(lc.near.h, lc.near.ppm, (k) => lc.near.e[k] % 2 !== 0), even = inT(lc.near.h, lc.near.ppm, (k) => lc.near.e[k] % 2 === 0);
    roF.set('in', `${fmt(-inT(lc.near.h, lc.near.ppm))} ppm`);
    roF.set('odd', Number.isFinite(odd) ? `${fmt(-odd)} ppm` : '–'); roF.set('even', Number.isFinite(even) ? `${fmt(-even)} ppm` : '–');
    const sec = inT(lc.sec.h, lc.sec.ppm);
    roF.set('sec', Number.isFinite(sec) ? `${fmt(-sec)} ppm` : '–');
  }
}
function maxWin(lc) { return Math.max(1, Math.max(...lc.near.h.map(Math.abs))); }

/* ───── check list ───── */
function checks(c) {
  const f = new Set((c.flags || '').split(',').filter(Boolean));
  const dfl = new Set((c.deep_flags || '').split(',').filter(Boolean));
  const out = [];
  const add = (name, state, val, expl) => out.push({ name, state, val, expl });
  const sc = (flag, name, val, expl) => add(name, f.has(flag) ? 'fail' : 'pass', val, expl);
  sc('low_sde', 'Periodogram peak stands out', `SDE ${fmt(c.sde, 1)} (≥ 9)`, 'How far the best period rises above all the other periods tried. Below 9, noise alone produces peaks this tall.');
  sc('low_snr', 'Transit signal-to-noise', `SNR ${fmt(c.snr, 1)} (≥ 7)`, 'Depth divided by its uncertainty, using every in-transit point.');
  sc('few_transits', 'At least 3 transits', `${c.n_transits} seen`, 'Two dips can line up by chance; three at equal spacing rarely do.');
  sc('odd_even', 'Odd and even transits match', `${fmt(c.odd_even_sigma, 1)}σ (≤ 3)`, 'Different depths for alternate dips point to an eclipsing binary at twice the period.');
  sc('secondary', 'Nothing at phase 0.5', `${fmt(c.secondary_sigma, 1)}σ (≤ 3)`, 'A dip half an orbit later is the companion being eclipsed: it glows, so it is a star.');
  sc('variable', 'Not a pulsating / spotted star', f.has('variable') ? 'brightening at phase 0.5' : 'ok', 'A brightening half an orbit later means a sine-like variation, not a transit.');
  sc('too_deep', 'Depth below 5%', `${fmt(c.depth_ppm)} ppm`, 'Planets almost never block more than a few percent of a Sun-like star.');
  sc('too_long', 'Duration fits the star', `${fmt(c.duration_h, 1)} h`, 'A transit cannot last much longer than a planet needs to cross the star at that orbital period.');
  sc('rp_too_large', 'Planet-sized', c.rp_earth != null ? `${fmt(c.rp_earth, 1)} R⊕ (≤ 25)` : 'no stellar radius', 'Depth × star size gives the size of the blocking body. Above ~2 Jupiters it is a star.');
  sc('common_period', 'Period not shared by many stars', f.has('common_period') ? 'shared' : 'unique', 'Spacecraft systematics (momentum dumps, scattered light) put the same fake period on many stars at once.');
  if (!c.deep_vetted) {
    add('Deep vetting', 'na', 'not run', c.status === 'known' || c.status === 'rejected' ? 'Deep checks run only on unmatched signals that pass every screen.' : 'Pending: runs after the sector finishes.');
    return out;
  }
  const cen = dfl.has('centroid_offset') || dfl.has('dip_off_target') ? 'fail' : dfl.has('centroid_shift_minor') || dfl.has('centroid_inconclusive') ? 'warn' : 'pass';
  add('Dip sits on this star (difference image)', cen, c.centroid_offset_px != null ? `${fmt(c.centroid_offset_px, 2)} px · ${fmt(c.centroid_offset_sigma, 1)}σ` : 'inconclusive',
    'We subtract in-transit pixels from out-of-transit pixels; the light that went missing should be centred on our star (one TESS pixel = 21″).');
  const nb = dfl.has('unresolved_neighbour') ? 'warn' : dfl.has('neighbour_could_mimic') ? 'warn' : 'pass';
  add('No neighbour bright enough to fake it', nb, c.gaia_n_can_mimic ? `${c.gaia_n_can_mimic} Gaia star(s)${c.gaia_mimic_min_sep != null ? `, closest ${fmt(c.gaia_mimic_min_sep, 1)}″` : ''}` : `none within 42″`,
    'A fainter star that dims completely could cause this depth only if it is bright enough. Ones closer than ~10″ can’t be ruled out with TESS pixels.');
  add('Not a catalogued eclipsing binary', dfl.has('known_eb') ? 'fail' : 'pass', c.eb_match || 'no match', 'Cross-match with the TESS (Prša et al. 2022) and Kepler eclipsing-binary catalogues, by TIC and by position.');
  add('Flat-bottomed, not V-shaped', dfl.has('v_shape') ? 'fail' : 'pass', c.v_ratio != null ? `V-ratio ${fmt(c.v_ratio, 2)} (box ≈ 1, V ≈ 0.4)` : '–', 'Depth in the middle third of the transit compared with the outer two thirds.');
  add('No secondary eclipse at any phase', dfl.has('secondary_scan') ? 'fail' : 'pass', c.sec_scan_sigma != null ? `best ${fmt(c.sec_scan_sigma, 1)}σ at phase ${fmt(c.sec_scan_phase, 2)} (≤ 5)` : '–', 'Eccentric binaries put the secondary away from phase 0.5, so we scan every phase.');
  add('Smaller than 2 Jupiters', dfl.has('rp_gt_2rjup') ? 'fail' : 'pass', c.rp_rjup != null ? `${fmt(c.rp_rjup, 2)} R_J` : 'no stellar radius', 'Using the TIC stellar radius.');
  return out;
}
function check(k) {
  const ic = { pass: '✓', warn: '!', fail: '✕', na: '–' }[k.state];
  const word = { pass: 'passed', warn: 'warning', fail: 'failed', na: 'not run' }[k.state];
  return el('div', { class: `check ${k.state}` }, el('span', { class: 'ic', 'aria-label': word, role: 'img' }, ic), el('b', {}, k.name), el('span', { class: 'val' }, k.val), el('p', {}, k.expl));
}

main().catch((e) => { document.getElementById('title').textContent = 'Could not load data'; document.getElementById('lede').textContent = String(e); });
