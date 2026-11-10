// Learn-page explorables. Each widget mounts into <div data-widget="…">, runs real (if simplified)
// physics in the browser, and comes with a "What to notice" list, in the spirit of the RF Field Notes explorers.
import { el, fmt, slider, segmented, toggle, readout, notice, surface, axes, palette, alpha, Loop, playButton, rng, starfield } from './ui.js';

starfield();
const R_EARTH_RSUN = 0.009168, AU_RSUN = 215.032;
const box = (s, l = 58, b = 40) => ({ x0: l, y0: 10, w: s.w - l - 12, h: s.h - 10 - b });
function frame(root, title, intro) {
  const viz = el('div', { class: 'viz' }), controls = el('aside', { class: 'controls', 'aria-label': `${title} controls` });
  root.append(el('h3', { class: 'wtitle' }, title), intro ? el('p', { class: 'kbd' }, intro) : '', el('div', { class: 'stage' }, viz, controls));
  return { viz, controls };
}
const canvasIn = (parent, label, aspect, minH = 200, maxH = 360, onR) => {
  const cv = el('canvas', { role: 'img', 'aria-label': label });
  parent.append(el('div', {}, cv));
  return surface(cv, { aspect, minH, maxH }, onR);
};

/* Quadratic limb darkening, roughly a Sun-like star in the TESS band. */
const LD = [0.4, 0.25];
const I = (r, ld) => { const mu = Math.sqrt(Math.max(0, 1 - r * r)); return ld ? 1 - LD[0] * (1 - mu) - LD[1] * (1 - mu) ** 2 : 1; };
const I_MEAN_LD = 1 - LD[0] / 3 - LD[1] / 6;
/** Fraction of starlight blocked by a dark disc of radius k (stellar radii) centred at distance z. */
function blocked(z, k, ld, n = 22) {
  if (z >= 1 + k) return 0;
  let s = 0;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    const dx = ((i + 0.5) / n * 2 - 1) * k, dy = ((j + 0.5) / n * 2 - 1) * k;
    if (dx * dx + dy * dy > k * k) continue;
    const r = Math.hypot(z + dx, dy);
    if (r < 1) s += I(r, ld);
  }
  const cell = (2 * k / n) ** 2;
  return (s * cell) / (Math.PI * (ld ? I_MEAN_LD : 1));
}

/* ───────── 1. Transit simulator ───────── */
function transitSim(root) {
  const { viz, controls } = frame(root, 'Transit simulator', 'Change the planet and the star and watch the dip respond. The light curve is computed live from the geometry.');
  const st = { rp: 2, rs: 1, p: 10, b: 0.3, ld: true, noise: false, tmag: 10 };
  const ro = readout([{ k: 'depth', label: 'Depth' }, { k: 'dur', label: 'Duration' }, { k: 'aR', label: 'a / R★' }, { k: 'prob', label: 'Chance aligned' },
    { k: 'snr', label: 'SNR, one sector' }, { k: 'n', label: 'Transits / sector' }]);
  let loop;
  const S1 = canvasIn(viz, 'Star and planet seen from Earth', 0.42, 200, 300, () => loop?.request());
  const S2 = canvasIn(viz, 'Light curve', 0.36, 200, 300, () => loop?.request());
  viz.append(ro.el);
  const geom = () => {
    const mass = Math.max(0.1, Math.min(2.5, st.rs));                          // main-sequence-ish M ≈ R
    const aAU = Math.cbrt(mass * (st.p / 365.25) ** 2);
    const aR = (aAU * AU_RSUN) / st.rs, k = (st.rp * R_EARTH_RSUN) / st.rs;
    const transits = st.b < 1 + k;
    const arg = Math.sqrt(Math.max(0, (1 + k) ** 2 - st.b ** 2)) / aR;
    const T = transits ? (st.p * 24 / Math.PI) * Math.asin(Math.min(1, arg)) : 0;   // hours
    return { aR, k, T, transits };
  };
  const curve = (g) => {
    const half = Math.max(g.T * 0.9, 1.5), N = 220, xs = [], ys = [];
    const vel = (2 * Math.PI * g.aR) / (st.p * 24);                            // stellar radii per hour
    for (let i = 0; i <= N; i++) {
      const t = -half + (2 * half * i) / N, x = vel * t;
      xs.push(t); ys.push(-blocked(Math.hypot(x, st.b), g.k, st.ld) * 1e6);
    }
    return { xs, ys, half, vel };
  };
  let cache = null;
  const update = () => {
    const g = geom(); cache = { g, c: curve(g) };
    const depth = -Math.min(...cache.c.ys);
    ro.set('depth', g.transits ? `${fmt(depth)} ppm` : 'misses the star', g.transits ? '' : 'badv');
    ro.set('dur', g.transits ? `${g.T.toFixed(2)} h` : '–');
    ro.set('aR', g.aR.toFixed(1)); ro.set('prob', `${(100 / g.aR).toFixed(2)}%`);
    const n = Math.floor(27 / st.p) || (st.p < 27 ? 1 : 0);
    ro.set('n', String(n), n < 3 ? 'warnv' : 'good');
    const sig1h = 60 * 10 ** (0.2 * (st.tmag - 10));                           // rough TESS noise per hour, ppm
    const snr = g.transits ? (depth / sig1h) * Math.sqrt(Math.max(g.T, 0.1) * n) : 0;
    ro.set('snr', g.transits ? snr.toFixed(1) : '–', snr >= 7 ? 'good' : 'badv');
    cache.sig = sig1h;
    loop?.request();
  };
  const r = rng(7);
  const noiseBank = Array.from({ length: 221 }, () => r.g());
  loop = new Loop((t) => {
    const P = palette();
    if (!cache) return;
    const { g, c } = cache;
    const phase = ((t * 0.25) % 1) * 2 - 1;                                    // -1..1 across the window
    const tNow = phase * c.half;
    { const { ctx, w, h } = S1;
      if (!w) return;
      ctx.clearRect(0, 0, w, h);
      const R = Math.min(h * 0.42, w * 0.22), cx = w / 2, cy = h / 2;
      const grd = ctx.createRadialGradient(cx, cy, 0, cx, cy, R);
      if (st.ld) { grd.addColorStop(0, '#fff3d6'); grd.addColorStop(0.7, '#ffcf87'); grd.addColorStop(1, '#c46a1c'); } else { grd.addColorStop(0, '#ffd59e'); grd.addColorStop(1, '#ffd59e'); }
      ctx.fillStyle = grd; ctx.beginPath(); ctx.arc(cx, cy, R, 0, 6.283); ctx.fill();
      ctx.strokeStyle = alpha(P.acc, 0.4); ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(0, cy - st.b * R); ctx.lineTo(w, cy - st.b * R); ctx.stroke(); ctx.setLineDash([]);
      const px = cx + c.vel * tNow * R, pr = Math.max(1.5, g.k * R);
      ctx.fillStyle = '#05070f'; ctx.strokeStyle = alpha(P.acc, 0.8); ctx.lineWidth = 1;
      ctx.beginPath(); ctx.arc(px, cy - st.b * R, pr, 0, 6.283); ctx.fill(); ctx.stroke();
      ctx.fillStyle = P.muted; ctx.font = '11px "JetBrains Mono", monospace'; ctx.textAlign = 'left';
      ctx.fillText(`b = ${st.b.toFixed(2)}`, 10, h - 10);
    }
    { const { ctx, w, h } = S2;
      ctx.clearRect(0, 0, w, h);
      const depth = Math.max(-Math.min(...c.ys), 50);
      const lo = st.noise ? -depth - 3 * cache.sig * 2.4 : -depth * 1.25, hi = st.noise ? 3 * cache.sig * 2.4 : depth * 0.25;
      const b_ = box(S2);
      const { X, Y } = axes(ctx, b_, { x: [-c.half, c.half], y: [lo, hi], xlabel: 'hours from mid-transit', ylabel: 'ppm', fy: (v) => fmt(v) });
      ctx.save(); ctx.beginPath(); ctx.rect(b_.x0, b_.y0, b_.w, b_.h); ctx.clip();
      if (st.noise) {                                                           // 10-min points, one transit
        const s10 = cache.sig * Math.sqrt(6);
        ctx.fillStyle = alpha(P.acc, 0.55);
        for (let i = 0; i < c.xs.length; i += 2) { const y = c.ys[i] + s10 * noiseBank[i]; ctx.beginPath(); ctx.arc(X(c.xs[i]), Y(y), 2, 0, 6.283); ctx.fill(); }
      }
      ctx.strokeStyle = P.acc2; ctx.lineWidth = 2.4; ctx.beginPath();
      c.xs.forEach((x, i) => (i ? ctx.lineTo(X(x), Y(c.ys[i])) : ctx.moveTo(X(x), Y(c.ys[i])))); ctx.stroke();
      ctx.strokeStyle = alpha(P.fg, 0.6); ctx.beginPath(); ctx.moveTo(X(tNow), b_.y0); ctx.lineTo(X(tNow), b_.y0 + b_.h); ctx.stroke();
      ctx.restore();
    }
  }, viz);
  const sRp = slider({ id: 'rp', label: 'Planet radius', min: 0.5, max: 25, log: true, value: st.rp, fmt: (v) => `${v.toFixed(1)} R⊕${v > 9 ? ` (${(v / 11.2).toFixed(2)} R_J)` : ''}`, onInput: (v) => { st.rp = v; update(); } });
  const sRs = slider({ id: 'rs', label: 'Star radius', min: 0.15, max: 3, log: true, value: st.rs, fmt: (v) => `${v.toFixed(2)} R☉`, onInput: (v) => { st.rs = v; update(); },
    hint: 'Same planet, smaller star: a much deeper dip. That is why small planets are easiest to find around M dwarfs.' });
  const sP = slider({ id: 'p', label: 'Orbital period', min: 0.5, max: 400, log: true, value: st.p, fmt: (v) => `${v.toFixed(2)} d`, onInput: (v) => { st.p = v; update(); } });
  const sB = slider({ id: 'b', label: 'Impact parameter b', min: 0, max: 1.2, step: 0.01, value: st.b, fmt: (v) => v.toFixed(2), onInput: (v) => { st.b = v; update(); },
    hint: '0 = across the middle of the star, 1 = grazing the edge.' });
  const tLD = toggle({ id: 'ld', label: 'Limb darkening (star dimmer at its edge)', value: true, onChange: (v) => { st.ld = v; update(); } });
  const tN = toggle({ id: 'noise', label: 'Add TESS-like noise (one transit)', value: false, onChange: (v) => { st.noise = v; update(); } });
  const sT = slider({ id: 'tmag', label: 'Star brightness (TESS mag)', min: 6, max: 16, step: 0.1, value: st.tmag, fmt: (v) => v.toFixed(1), onInput: (v) => { st.tmag = v; update(); },
    hint: 'Fainter stars are noisier: 5 magnitudes fainter is 10× the noise.' });
  const presets = segmented({ id: 'preset', label: 'Try', value: 'custom', options: [{ value: 'earth', label: 'Earth/Sun' }, { value: 'hotj', label: 'Hot Jupiter' }, { value: 'mdw', label: 'Super-Earth, M dwarf' }, { value: 'custom', label: 'Custom' }],
    onChange: (v) => {
      const p = { earth: [1, 1, 365.25, 0.2], hotj: [12, 1.1, 3.5, 0.3], mdw: [1.6, 0.25, 4, 0.2] }[v];
      if (!p) return;
      [st.rp, st.rs, st.p, st.b] = p; sRp.set(st.rp); sRs.set(st.rs); sP.set(st.p); sB.set(st.b); update();
    } });
  controls.append(el('div', {}, playButton(loop)), presets.el, el('h3', {}, 'Planet'), sRp.el, sP.el, sB.el, el('h3', {}, 'Star'), sRs.el, tLD.el, el('h3', {}, 'Data'), tN.el, sT.el);
  root.append(notice([
    '<b>Depth ≈ (R<sub>planet</sub>/R<sub>star</sub>)²</b>: an Earth in front of the Sun is ~84 ppm, a Jupiter ~1%.',
    '<b>Raise b toward 1</b> and the transit gets shorter and shallower, and its bottom rounds into a V. Past 1 + R<sub>p</sub>/R★ there is no transit at all.',
    '<b>Limb darkening</b> rounds the bottom of the dip: the planet blocks brighter light near the centre of the star than at the edge.',
    '<b>Long periods mean few transits</b>: TESS watches most stars for 27 days, so anything slower than ~9 days shows fewer than 3 dips, and our search won’t trust it.',
    '<b>Chance aligned</b> is R★/a: for an Earth twin it is 0.5%, which is why surveys have to watch so many stars.',
  ]));
  update();
}

/* ───────── 2. BLS by hand ───────── */
function blsDemo(root) {
  const { viz, controls } = frame(root, 'Find the period yourself', 'A fake 27-day light curve with a planet buried in noise. Each transit is invisible on its own. Drag the trial period until the dips stack up.');
  const r = rng(42), TRUE_P = 3.4172, T0 = 1.3, DUR = 2.6 / 24, DEPTH = 700, SIG = 1100;
  const t = [], y = [];
  for (let x = 0; x < 27; x += 10 / 1440) {
    if (x > 13.1 && x < 14.3) continue;                                         // data downlink gap
    const ph = ((x - T0) % TRUE_P + TRUE_P) % TRUE_P, d = Math.min(ph, TRUE_P - ph);
    t.push(x); y.push((d < DUR / 2 ? -DEPTH : 0) + SIG * r.g() + 300 * Math.sin(x / 4.1));
  }
  // remove the slow wiggle the way the pipeline would (running median, ~0.75 d)
  const flat = y.map((_, i) => { const w = []; for (let j = Math.max(0, i - 54); j < Math.min(y.length, i + 54); j++) w.push(y[j]); w.sort((a, b) => a - b); return y[i] - w[w.length >> 1]; });
  const NB = (P) => Math.max(8, Math.round(P / DUR));
  const fold = (P) => {
    const nb = NB(P), sum = new Float64Array(nb), cnt = new Uint32Array(nb);
    for (let i = 0; i < t.length; i++) { const ph = (((t[i] - T0 + DUR / 2) / P) % 1 + 1) % 1; const j = Math.floor(ph * nb); sum[j] += flat[i]; cnt[j]++; }
    let best = 0, bj = 0;
    for (let j = 0; j < nb; j++) if (cnt[j]) { const z = -(sum[j] / cnt[j]) / (SIG / Math.sqrt(cnt[j])); if (z > best) { best = z; bj = j; } }
    return { best, bj, nb, sum, cnt };
  };
  const grid = [], power = [];
  for (let P = 1; P <= 9; P *= 1.0004) { grid.push(P); power.push(fold(P).best); }
  const st = { p: 2.3 };
  const S0 = canvasIn(viz, 'Unfolded light curve', 0.25, 150, 220, () => draw());
  const S1 = canvasIn(viz, 'Folded at the trial period', 0.34, 200, 300, () => draw());
  const S2 = canvasIn(viz, 'Periodogram: dip strength at every trial period', 0.26, 160, 240, () => draw());
  const ro = readout([{ k: 'p', label: 'Trial period' }, { k: 'z', label: 'Best dip' }, { k: 'msg', label: 'Status' }]);
  viz.append(ro.el);
  const P_ = palette();
  function draw() {
    const P = P_;
    { const { ctx, w, h } = S0; if (!w) return; ctx.clearRect(0, 0, w, h);
      const b = box(S0, 58, 34); const { X, Y } = axes(ctx, b, { x: [0, 27], y: [-4500, 4500], xlabel: 'days', ylabel: 'ppm', ny: 2, fy: (v) => fmt(v) });
      ctx.fillStyle = alpha(P.acc, 0.35); for (let i = 0; i < t.length; i += 2) ctx.fillRect(X(t[i]), Y(flat[i]), 1.5, 1.5);
      ctx.strokeStyle = alpha(P.acc2, 0.7);
      for (let x = T0; x < 27; x += st.p) { ctx.beginPath(); ctx.moveTo(X(x), b.y0); ctx.lineTo(X(x), b.y0 + 10); ctx.stroke(); }
    }
    const f = fold(st.p);
    { const { ctx, w, h } = S1; ctx.clearRect(0, 0, w, h);
      const b = box(S1); const { X, Y } = axes(ctx, b, { x: [0, 1], y: [-1800, 900], xlabel: 'phase (fraction of the trial period)', ylabel: 'ppm', fy: (v) => fmt(v) });
      ctx.save(); ctx.beginPath(); ctx.rect(b.x0, b.y0, b.w, b.h); ctx.clip();
      ctx.fillStyle = alpha(P.acc, 0.18);
      for (let i = 0; i < t.length; i++) { const ph = (((t[i] - T0 + DUR / 2) / st.p) % 1 + 1) % 1; ctx.fillRect(X(ph), Y(flat[i]), 1.5, 1.5); }
      ctx.strokeStyle = P.acc2; ctx.lineWidth = 2.4; ctx.beginPath();
      for (let j = 0; j < f.nb; j++) { if (!f.cnt[j]) continue; const xx = X((j + 0.5) / f.nb), yy = Y(f.sum[j] / f.cnt[j]); j ? ctx.lineTo(xx, yy) : ctx.moveTo(xx, yy); }
      ctx.stroke(); ctx.restore();
    }
    { const { ctx, w, h } = S2; ctx.clearRect(0, 0, w, h);
      const b = box(S2, 58, 34); const lx = (p) => Math.log(p);
      const { X, Y } = axes(ctx, b, { x: [lx(1), lx(9)], y: [0, Math.max(...power) * 1.1], xlabel: 'trial period (days, log)', ylabel: 'σ', ny: 3, fx: (v) => Math.exp(v).toFixed(1) });
      ctx.strokeStyle = P.sub; ctx.lineWidth = 1; ctx.beginPath(); grid.forEach((p, i) => (i ? ctx.lineTo(X(lx(p)), Y(power[i])) : ctx.moveTo(X(lx(p)), Y(power[i])))); ctx.stroke();
      ctx.strokeStyle = P.acc2; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(X(lx(st.p)), b.y0); ctx.lineTo(X(lx(st.p)), b.y0 + b.h); ctx.stroke();
      S2.box = b; S2.X = X;
    }
    const rel = st.p / TRUE_P, near = (m) => Math.abs(rel - m) < 0.0015;
    ro.set('p', `${st.p.toFixed(4)} d`); ro.set('z', `${f.best.toFixed(1)} σ`, f.best > 7 ? 'good' : '');
    ro.set('msg', near(1) ? 'Locked on: the true period!' : near(0.5) ? 'Half the period: every other “transit” is empty' : near(2) ? 'Double: only every 2nd transit stacks'
      : near(1 / 3) || near(3) ? 'A harmonic: close, not it' : f.best > 5 ? 'Something is lining up…' : 'Noise. Keep sliding.', near(1) ? 'good' : near(0.5) || near(2) ? 'warnv' : '');
  }
  const sP = slider({ id: 'tp', label: 'Trial period', min: 1, max: 9, log: true, value: st.p, fmt: (v) => `${v.toFixed(3)} d`, onInput: (v) => { st.p = v; draw(); } });
  const sF = slider({ id: 'fine', label: 'Fine tune', min: -0.02, max: 0.02, step: 0.0005, value: 0, fmt: (v) => `${v >= 0 ? '+' : ''}${v.toFixed(4)} d`,
    onInput: (v) => { st.p = Math.max(1, Math.min(9, sP.get() + v)); draw(); }, hint: 'Near the right period, tiny changes matter: an error of 0.01 d smears 8 transits by 2 hours.' });
  const best = el('button', { type: 'button', class: 'btn' }, 'Jump to the tallest peak (what BLS does)');
  best.addEventListener('click', () => { const i = power.indexOf(Math.max(...power)); st.p = grid[i]; sP.set(st.p); sF.set(0); draw(); });
  S2.canvas.addEventListener('click', (e) => {
    const rc = S2.canvas.getBoundingClientRect(), b = S2.box; if (!b) return;
    const x = Math.min(1, Math.max(0, (e.clientX - rc.left - b.x0) / b.w));
    st.p = Math.exp(Math.log(1) + x * (Math.log(9) - Math.log(1))); sP.set(st.p); sF.set(0); draw();
  });
  controls.append(sP.el, sF.el, best, el('small', { class: 'kbd' }, 'Tip: click the periodogram to jump there.'));
  root.append(notice([
    '<b>Folding</b> means cutting the light curve into pieces one trial-period long and stacking them. Only at the right period do the dips land on top of each other.',
    '<b>Noise shrinks as √N</b>: eight stacked transits are √8 ≈ 2.8× clearer than one.',
    '<b>Aliases</b>: at half the true period, half the stacked “transits” are empty, so the dip is diluted; at double, you only stack every other one. Both show up as smaller peaks.',
    '<b>BLS</b> does this for tens of thousands of trial periods and durations per star, then asks how far the best peak stands above the rest (the SDE).',
  ]));
  draw();
}

/* ───────── 3. Planet or eclipsing binary? ───────── */
function overlap(d, r1, r2) {                                                   // area of two overlapping circles
  if (d >= r1 + r2) return 0;
  if (d <= Math.abs(r1 - r2)) return Math.PI * Math.min(r1, r2) ** 2;
  const a = r1 * r1 * Math.acos((d * d + r1 * r1 - r2 * r2) / (2 * d * r1)), b = r2 * r2 * Math.acos((d * d + r2 * r2 - r1 * r1) / (2 * d * r2));
  return a + b - 0.5 * Math.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2));
}
function ebDemo(root) {
  const { viz, controls } = frame(root, 'Planet or eclipsing binary?', 'The commonest impostor is a pair of stars eclipsing each other. Pick a scenario and look for the tells our vetting uses.');
  const SC = {
    planet: { r2: 0.1, s2: 0, b: 0.3, blend: 0, label: 'Planet' },
    eb: { r2: 0.8, s2: 0.55, b: 0.2, blend: 0, label: 'Eclipsing binary' },
    graze: { r2: 0.7, s2: 0.45, b: 1.45, blend: 0, label: 'Grazing binary' },
    bg: { r2: 0.9, s2: 0.8, b: 0.5, blend: 0.985, label: 'Background binary (blend)' },
  };
  const st = { ...SC.eb, half: false };
  const S1 = canvasIn(viz, 'The system, to scale', 0.3, 150, 220, () => draw());
  const S2 = canvasIn(viz, 'One full orbit', 0.32, 180, 280, () => draw());
  const S3 = canvasIn(viz, 'Zoom on the primary eclipse', 0.32, 180, 280, () => draw());
  const ro = readout([{ k: 'd1', label: 'Primary depth' }, { k: 'd2', label: 'Secondary depth' }, { k: 'v', label: 'V-ratio' }, { k: 'oe', label: 'Odd vs even' }, { k: 'verdict', label: 'Our checks say' }]);
  viz.append(ro.el);
  const model = (phase) => {                                                    // phase 0 = primary eclipse
    const ang = 2 * Math.PI * phase, aR = 4;                                   // separation in primary radii
    const x = aR * Math.sin(ang), front = Math.cos(ang) > 0;                   // star 2 in front near phase 0
    const d = Math.hypot(x, st.b);
    const L1 = Math.PI, L2 = Math.PI * st.r2 ** 2 * st.s2, tot = L1 + L2;
    const ov = overlap(d, 1, st.r2);
    let lost = front ? ov : ov * st.s2;                                        // star 1 surface brightness = 1
    const frac = lost / tot;
    return -frac * (1 - st.blend) * 1e6;
  };
  const P = palette();
  function draw() {
    { const { ctx, w, h } = S1; if (!w) return; ctx.clearRect(0, 0, w, h);
      const R = Math.min(h * 0.3, w * 0.12), cy = h / 2;
      const star = (x, r, s) => { ctx.fillStyle = s > 0 ? `rgba(255,${180 + 50 * s},${100 + 80 * s},1)` : '#05070f'; ctx.beginPath(); ctx.arc(x, cy, r, 0, 6.283); ctx.fill(); if (s <= 0) { ctx.strokeStyle = alpha(P.acc, 0.8); ctx.stroke(); } };
      if (st.blend) {
        ctx.fillStyle = '#fff3d6'; ctx.beginPath(); ctx.arc(w * 0.28, cy, R * 1.3, 0, 6.283); ctx.fill();
        ctx.fillStyle = P.muted; ctx.font = '11px "JetBrains Mono", monospace'; ctx.textAlign = 'center';
        ctx.fillText('our target', w * 0.28, cy + R * 1.3 + 14); ctx.fillText('faint pair behind it', w * 0.7, cy + R * 0.9 + 14);
        star(w * 0.66, R * 0.45, 1); star(w * 0.66 + R * 0.45 * 1.6, R * 0.45 * st.r2, st.s2);
      } else { star(w * 0.45, R, 1); star(w * 0.45 + R * 2.2, R * Math.max(st.r2, 0.06), st.s2); }
    }
    const ph = [], fl = [];
    for (let i = 0; i <= 600; i++) { ph.push(i / 600 - 0.25); fl.push(model(i / 600 - 0.25)); }
    const d1 = -Math.min(...fl.slice(100, 200).concat(fl.slice(0, 1))), d2 = -Math.min(...fl.slice(400, 500));
    const lo = -Math.max(d1, 200) * 1.25;
    { const { ctx, w, h } = S2; ctx.clearRect(0, 0, w, h);
      const b = box(S2); const { X, Y } = axes(ctx, b, { x: [-0.25, 0.75], y: [lo, -lo * 0.15], xlabel: 'orbital phase', ylabel: 'ppm', fy: (v) => fmt(v) });
      ctx.strokeStyle = P.acc2; ctx.lineWidth = 2.2; ctx.beginPath(); ph.forEach((p, i) => (i ? ctx.lineTo(X(p), Y(fl[i])) : ctx.moveTo(X(p), Y(fl[i])))); ctx.stroke();
      if (st.half) { ctx.fillStyle = P.muted; ctx.font = '11px "JetBrains Mono", monospace'; ctx.fillText('BLS folds both eclipses together →', X(0.12), Y(lo * 0.6)); }
    }
    let ecl = 0; for (let p = 0; p < 0.25; p += 0.0005) if (model(p) < -1e-3 || model(p + 0.5) < -1e-3) ecl = p;
    { const { ctx, w, h } = S3; ctx.clearRect(0, 0, w, h);
      const b = box(S3); const z = Math.max(0.01, Math.min(0.25, ecl * 1.6));
      const { X, Y } = axes(ctx, b, { x: [-z, z], y: [lo, -lo * 0.15], xlabel: 'phase around the eclipse', ylabel: 'ppm', fy: (v) => fmt(v), fx: (v) => v.toFixed(3) });
      ctx.lineWidth = 2.2;
      const series = st.half ? [[0, P.acc2, 'odd'], [0.5, P.sub, 'even']] : [[0, P.acc2, '']];
      for (const [off, col] of series) { ctx.strokeStyle = col; ctx.beginPath(); for (let i = 0; i <= 200; i++) { const p = -z + (2 * z * i) / 200; const yv = model(p + off); i ? ctx.lineTo(X(p), Y(yv)) : ctx.moveTo(X(p), Y(yv)); } ctx.stroke(); }
    }
    // V-ratio: mean depth of the middle third / outer two thirds of the eclipse (as in our deep vetting)
    let half = 0; for (let p = 0; p < 0.25; p += 0.0005) if (model(p) < -1e-3) half = p;
    const avg = (a, b_) => { let s = 0, n = 0; for (let p = a; p <= b_; p += half / 60) { s += model(p); n++; } return -s / n; };
    const v = half ? avg(0.5 * half, half) / avg(0, half / 3) : NaN;
    ro.set('d1', `${fmt(d1)} ppm`); ro.set('d2', d2 > 1 ? `${fmt(d2)} ppm` : 'none', d2 > 0.1 * d1 && d2 > 50 ? 'badv' : 'good');
    ro.set('v', Number.isFinite(v) ? v.toFixed(2) : '–', v < 0.45 ? 'badv' : 'good');
    const oe = st.half ? Math.abs(d1 - d2) / Math.max(d1, 1) : 0;
    ro.set('oe', st.half ? `${fmt(d1)} vs ${fmt(d2)} ppm` : 'fold at the true period', st.half && oe > 0.15 ? 'badv' : '');
    const tells = [];
    if (d2 > 0.1 * d1 && d2 > 50) tells.push(st.half ? 'odd ≠ even' : 'secondary eclipse');
    if (v < 0.45) tells.push('V-shaped');
    if (st.blend) tells.push('difference image off target');
    ro.set('verdict', tells.length ? tells.join(' · ') : 'consistent with a planet', tells.length ? 'badv' : 'good');
  }
  const seg = segmented({ id: 'scn', label: 'Scenario', value: 'eb', options: Object.entries(SC).map(([k, v]) => ({ value: k, label: v.label })),
    onChange: (k) => { Object.assign(st, SC[k]); sR.set(st.r2); sS.set(st.s2); sB.set(st.b); sL.set(st.blend * 100); draw(); } });
  const sR = slider({ id: 'r2', label: 'Companion radius', min: 0.03, max: 1, step: 0.01, value: st.r2, fmt: (v) => `${v.toFixed(2)} R★`, onInput: (v) => { st.r2 = v; draw(); } });
  const sS = slider({ id: 's2', label: 'Companion surface brightness', min: 0, max: 1, step: 0.01, value: st.s2, fmt: (v) => (v ? v.toFixed(2) : 'dark (planet)'), onInput: (v) => { st.s2 = v; draw(); },
    hint: 'A planet is dark, so it never makes a secondary eclipse TESS can see. A star glows.' });
  const sB = slider({ id: 'eb_b', label: 'Impact parameter', min: 0, max: 1.9, step: 0.01, value: st.b, fmt: (v) => v.toFixed(2), onInput: (v) => { st.b = v; draw(); } });
  const sL = slider({ id: 'blend', label: 'Light from other stars', min: 0, max: 99.9, step: 0.1, value: st.blend * 100, fmt: (v) => `${v.toFixed(1)}%`, onInput: (v) => { st.blend = v / 100; draw(); },
    hint: 'A bright foreground star dilutes a deep eclipse into a planet-sized dip.' });
  const tH = toggle({ id: 'half', label: 'Fold at half the true period (what BLS often finds)', value: false, onChange: (v) => { st.half = v; draw(); } });
  controls.append(seg.el, sR.el, sS.el, sB.el, sL.el, tH.el);
  root.append(notice([
    '<b>Secondary eclipse</b>: when the companion glows, it is eclipsed half an orbit later. Planets are too faint for that to show.',
    '<b>Equal stars give equal eclipses</b>, so BLS often locks on half the true period. Then odd and even “transits” have different depths if the stars differ at all.',
    '<b>Grazing</b> eclipses never have a flat bottom: the V-ratio drops well below 1.',
    '<b>Blends</b> can look perfect in the light curve. Only the pixels give them away: next section.',
  ]));
  draw();
}

/* ───────── 4. Whose dip is it? (difference-image centroids) ───────── */
function centroidDemo(root) {
  const { viz, controls } = frame(root, 'Whose dip is it?', 'An 11×11 TESS pixel stamp around our target, with one neighbour. Choose which star actually dims and see what the difference image says.');
  const st = { sep: 30, dm: 2, who: 'target', depth: 2000, psf: 0.9 };
  const N = 11, C0 = 5;
  const S = canvasIn(viz, 'Pixel images', 0.5, 220, 360, () => draw());
  const ro = readout([{ k: 'need', label: 'Neighbour must dim by' }, { k: 'off', label: 'Centroid offset' }, { k: 'verdict', label: 'Verdict' }]);
  viz.append(ro.el);
  const P = palette();
  const img = (sc1, sc2) => {
    const sepPx = st.sep / 21, nx = C0 + sepPx * Math.cos(0.6), ny = C0 + sepPx * Math.sin(0.6), f2 = 10 ** (-0.4 * st.dm);
    const a = [];
    for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
      const g = (cx, cy) => Math.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * st.psf ** 2));
      a.push(sc1 * g(C0, C0) + sc2 * f2 * g(nx, ny));
    }
    return { a, nx, ny, f2 };
  };
  const centroid = (a) => {                                                     // 3×3 around the brightest pixel, as in our pipeline
    let m = -Infinity, mi = 0; a.forEach((v, i) => { if (v > m) { m = v; mi = i; } });
    const iy = Math.floor(mi / N), ix = mi % N; let s = 0, sx = 0, sy = 0;
    for (let y = Math.max(0, iy - 1); y <= Math.min(N - 1, iy + 1); y++) for (let x = Math.max(0, ix - 1); x <= Math.min(N - 1, ix + 1); x++) { const v = Math.max(0, a[y * N + x]); s += v; sx += v * x; sy += v * y; }
    return { x: sx / s, y: sy / s, ix, iy };
  };
  function draw() {
    const { ctx, w, h } = S; if (!w) return; ctx.clearRect(0, 0, w, h);
    const out = img(1, 1);
    const tot = 1 + out.f2, lost = (st.depth / 1e6) * tot;                     // flux the dip removes
    const need = lost / out.f2;                                                 // fraction of the neighbour that must vanish
    const possible = st.who === 'target' || need <= 1;
    const dimT = st.who === 'target' ? lost : 0, dimN = st.who === 'target' ? 0 : Math.min(1, need) * out.f2;
    const diff = img(dimT, dimN / out.f2).a;
    const pane = Math.min((w - 30) / 2, h - 34), cell = pane / N;
    const draw1 = (a, x0, title, cen, mark) => {
      const mx = Math.max(...a) || 1;
      for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
        const v = Math.max(0, a[y * N + x] / mx) ** 0.6;
        ctx.fillStyle = `rgb(${Math.round(30 + 225 * v)},${Math.round(20 + 150 * v * v)},${Math.round(60 + 40 * v)})`;
        ctx.fillRect(x0 + x * cell, 22 + y * cell, cell - 1, cell - 1);
      }
      ctx.fillStyle = P.fg; ctx.font = '12px Inter, system-ui'; ctx.textAlign = 'left'; ctx.fillText(title, x0, 14);
      const px = (x) => x0 + (x + 0.5) * cell, py = (y) => 22 + (y + 0.5) * cell;
      ctx.lineWidth = 2;
      ctx.strokeStyle = P.acc; ctx.beginPath(); ctx.moveTo(px(C0) - 8, py(C0)); ctx.lineTo(px(C0) + 8, py(C0)); ctx.moveTo(px(C0), py(C0) - 8); ctx.lineTo(px(C0), py(C0) + 8); ctx.stroke();
      ctx.strokeStyle = alpha(P.fg, 0.6); ctx.beginPath(); ctx.arc(px(out.nx), py(out.ny), 5, 0, 6.283); ctx.stroke();
      if (mark) { ctx.strokeStyle = P.acc2; ctx.beginPath(); ctx.moveTo(px(cen.x) - 6, py(cen.y) - 6); ctx.lineTo(px(cen.x) + 6, py(cen.y) + 6); ctx.moveTo(px(cen.x) + 6, py(cen.y) - 6); ctx.lineTo(px(cen.x) - 6, py(cen.y) + 6); ctx.stroke(); }
    };
    const cOut = centroid(out.a), cDiff = centroid(diff);
    draw1(out.a, 10, 'Out of transit', cOut, false);
    if (possible) draw1(diff, 20 + pane, 'Light lost in transit', cDiff, true);
    else { ctx.fillStyle = P.bad; ctx.font = '13px Inter, system-ui'; ctx.fillText('Impossible: the neighbour is too faint', 20 + pane, 40); ctx.fillText('to cause this dip even if it vanished.', 20 + pane, 60); }
    const off = Math.hypot(cDiff.x - cOut.x, cDiff.y - cOut.y), peak = Math.hypot(cDiff.ix - C0, cDiff.iy - C0);
    ro.set('need', need > 1 ? `${(need * 100).toFixed(0)}% (impossible)` : `${(need * 100).toFixed(need < 0.01 ? 2 : 1)}%`, need > 1 ? 'good' : need > 0.5 ? 'warnv' : 'badv');
    ro.set('off', possible ? `${off.toFixed(2)} px = ${(off * 21).toFixed(0)}″` : '–');
    const flagged = possible && (off > 1 || peak > 2);
    const tooClose = st.sep < 10 && need <= 1;
    ro.set('verdict', !possible ? 'neighbour ruled out' : flagged ? 'flag: dip off target' : tooClose ? 'can’t tell: neighbour unresolved' : off > 0.5 ? 'warning: small shift' : 'on target',
      !possible || (!flagged && !tooClose && off <= 0.5) ? 'good' : flagged ? 'badv' : 'warnv');
  }
  const seg = segmented({ id: 'who', label: 'Which star really dims?', value: st.who, options: [{ value: 'target', label: 'Our target' }, { value: 'nb', label: 'The neighbour' }], onChange: (v) => { st.who = v; draw(); } });
  const sS = slider({ id: 'sep', label: 'Neighbour separation', min: 2, max: 90, step: 1, value: st.sep, fmt: (v) => `${v}″ (${(v / 21).toFixed(1)} px)`, onInput: (v) => { st.sep = v; draw(); } });
  const sM = slider({ id: 'dm', label: 'Neighbour is fainter by', min: 0, max: 8, step: 0.1, value: st.dm, fmt: (v) => `${v.toFixed(1)} mag (${(100 * 10 ** (-0.4 * v)).toFixed(v > 5 ? 2 : 1)}% as bright)`, onInput: (v) => { st.dm = v; draw(); } });
  const sD = slider({ id: 'cdepth', label: 'Observed depth', min: 100, max: 20000, log: true, value: st.depth, fmt: (v) => `${fmt(v)} ppm`, onInput: (v) => { st.depth = v; draw(); } });
  controls.append(seg.el, sS.el, sM.el, sD.el);
  root.append(notice([
    '<b>One TESS pixel is 21 arcseconds</b>, so several stars often share the same few pixels and their light is added together.',
    '<b>The difference image</b> (before/after minus during) shows where light went missing. If it is centred on the neighbour, the dip is the neighbour’s.',
    '<b>Brightness sets the limit</b>: a neighbour fainter by Δm can fake a depth of at most 10<sup>−Δm/2.5</sup>. Make it faint enough and it is ruled out without any pixels.',
    '<b>Closer than about 10″</b> the two blur together and even the difference image can’t separate them. Those we only warn about; telescopes with sharper images settle it.',
  ]));
  draw();
}

/* ───────── 5. Quiz ───────── */
const QUIZ = [
  { q: 'A planet half the radius of another, around the same star, makes a transit that is…', a: ['half as deep', 'a quarter as deep', 'twice as long', 'the same depth'], c: 1, why: 'Depth scales with the blocked area, (Rp/R★)². Halve the radius, quarter the depth.' },
  { q: 'Why do we fold a light curve at a trial period?', a: ['To remove the star’s rotation', 'To stack many faint transits so they add up', 'To convert time to frequency', 'To find the star’s radius'], c: 1, why: 'Stacking N transits beats the noise down by √N; only the right period lines them up.' },
  { q: 'A dip shows up halfway between transits, at phase 0.5. The most likely explanation is…', a: ['a second planet', 'an eclipsing binary star', 'limb darkening', 'a data gap'], c: 1, why: 'That is the secondary eclipse: the companion is being hidden, so it gives off light, so it is a star.' },
  { q: 'Every other transit is noticeably deeper. What happened?', a: ['The planet is spinning', 'Two different stars eclipsing in turn, folded at half the true period', 'Star spots', 'TESS changed cameras'], c: 1, why: 'Odd/even depth differences are a classic binary tell: BLS stacked primary and secondary eclipses together.' },
  { q: 'The difference image is centred one and a half pixels from our star. What do we conclude?', a: ['The planet is large', 'The dip probably comes from a neighbouring star', 'The period is wrong', 'Nothing, pixels don’t matter'], c: 1, why: 'Light was lost from somewhere else. One pixel is 21″, enough to hold other stars.' },
  { q: 'A neighbour is 5 magnitudes fainter than our target (1% as bright). Can it fake a 2% dip?', a: ['Yes, easily', 'Only if it is a giant', 'No: even vanishing completely it removes at most 1%', 'Only at phase 0.5'], c: 2, why: 'A star can remove at most all of its own light. 1% < 2%, so it is ruled out.' },
  { q: 'A transit has a sharp V shape instead of a flat bottom. That suggests…', a: ['a small rocky planet', 'a grazing eclipse, often by a star', 'a very long period', 'a bright host star'], c: 1, why: 'Once a small body is fully on the disc the light curve flattens. A V means it never fully gets on: a graze, typical of binaries.' },
  { q: 'Our search re-found 48% of the known TOIs in sector 70. What is that number for?', a: ['It is the discovery rate', 'It measures our sensitivity, so we know what we miss', 'It is the false-positive rate', 'It is the fraction of stars with planets'], c: 1, why: 'Recovering known planets blind is how a search measures its own completeness.' },
  { q: 'Something passes every check on this site. Is it a planet?', a: ['Yes', 'Only if the SDE is above 20', 'Not yet: it needs human vetting, more data and follow-up observations', 'Yes, if it is in the TESS input catalogue'], c: 2, why: 'Automated checks only remove the obvious impostors. Confirmation needs independent observations, which is why nothing here is called a discovery.' },
];
function quiz(root) {
  const KEY = 'transit-hunt-quiz-v1';
  let saved = {}; try { saved = JSON.parse(localStorage.getItem(KEY) || '{}'); } catch { /* private mode */ }
  const score = el('p', { class: 'kbd', role: 'status', 'aria-live': 'polite' });
  const list = el('ol', { class: 'quiz' });
  const sync = () => { const n = Object.values(saved).filter((v) => v.ok).length, a = Object.keys(saved).length; score.textContent = `${n} right of ${a} answered · ${QUIZ.length} questions`; };
  QUIZ.forEach((item, qi) => {
    const fb = el('p', { class: 'fb', 'aria-live': 'polite' });
    const opts = el('div', { class: 'opts', role: 'group', 'aria-label': `Question ${qi + 1} answers` });
    const btns = item.a.map((txt, ai) => {
      const b = el('button', { type: 'button', class: 'opt' }, txt);
      b.addEventListener('click', () => { saved[qi] = { pick: ai, ok: ai === item.c }; try { localStorage.setItem(KEY, JSON.stringify(saved)); } catch { /* ignore */ } paint(); sync(); });
      opts.append(b); return b;
    });
    const paint = () => {
      const s = saved[qi];
      btns.forEach((b, ai) => { b.classList.toggle('right', !!s && ai === item.c); b.classList.toggle('wrong', !!s && ai === s.pick && !s.ok); b.setAttribute('aria-pressed', String(!!s && ai === s.pick)); });
      fb.innerHTML = s ? `<b>${s.ok ? 'Right.' : 'Not quite.'}</b> ${item.why}` : '';
    };
    paint();
    list.append(el('li', {}, el('p', { class: 'q' }, item.q), opts, fb));
  });
  const reset = el('button', { type: 'button', class: 'btn' }, 'Start over');
  reset.addEventListener('click', () => { saved = {}; try { localStorage.removeItem(KEY); } catch { /* ignore */ } list.replaceChildren(); root.replaceChildren(); quiz(root); });
  root.append(score, list, reset);
  sync();
}

const W = { transit: transitSim, bls: blsDemo, eb: ebDemo, centroid: centroidDemo, quiz };
for (const node of document.querySelectorAll('[data-widget]')) {
  try { W[node.dataset.widget]?.(node); } catch (e) { node.append(el('p', { class: 'kbd' }, 'This interactive could not start in your browser.')); console.error(e); }
}
