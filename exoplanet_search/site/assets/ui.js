// Small shared toolkit for the interactive pieces of the site. Plain ES module, no build step.
// Patterns (slider with live <output>, segmented radio group, readout grid, HiDPI canvas that
// tracks its container, an animation loop that respects reduced motion and pauses off-screen)
// follow the "Field Notes" RF explorers.

export const $ = (s, r = document) => r.querySelector(s);
export const $$ = (s, r = document) => [...r.querySelectorAll(s)];
export const reducedMotion = () => matchMedia('(prefers-reduced-motion: reduce)').matches;

/** el('div', {class: 'x', onclick: f}, child, 'text', ...) */
export function el(tag, attrs = {}, ...kids) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v == null || v === false) continue;
    if (k === 'class') e.className = v;
    else if (k === 'html') e.innerHTML = v;
    else if (k === 'text') e.textContent = v;
    else if (k.startsWith('on') && typeof v === 'function') e.addEventListener(k.slice(2), v);
    else e.setAttribute(k, v === true ? '' : v);
  }
  for (const c of kids.flat()) if (c != null && c !== false) e.append(c instanceof Node ? c : document.createTextNode(String(c)));
  return e;
}

export const fmt = (x, d = 0) => (x == null || Number.isNaN(x) ? '–'
  : Number(x).toLocaleString('en-US', { maximumFractionDigits: d, minimumFractionDigits: d }));

const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(`--${n}`).trim();
let pal = null;
export function palette() {
  if (!pal) pal = { fg: css('fg'), muted: css('muted'), acc: css('acc'), acc2: css('acc2'), ok: css('ok'), warn: css('warn'),
    bad: css('bad'), sub: css('sub'), line: 'rgba(124,196,255,.14)', grid: 'rgba(124,196,255,.08)', bg: css('bg2') };
  return pal;
}
export function alpha(hex, a) {
  const h = hex.replace('#', '');
  const n = parseInt(h.length === 3 ? [...h].map((x) => x + x).join('') : h, 16);
  return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`;
}

let uid = 0;
/** Range slider with a live value, aria-valuetext and an optional log scale. */
export function slider(o) {
  const id = `in-${o.id}-${++uid}`, N = 1000;
  const toPos = (v) => (o.log ? (Math.log(v / o.min) / Math.log(o.max / o.min)) * N : v);
  const fromPos = (p) => (o.log ? o.min * (o.max / o.min) ** (p / N) : p);
  const input = el('input', { type: 'range', id, min: o.log ? 0 : o.min, max: o.log ? N : o.max,
    step: o.log ? 1 : (o.step ?? 'any'), value: toPos(o.value), 'aria-describedby': o.hint ? `${id}-h` : null });
  const out = el('output', { for: id });
  const f = o.fmt ?? ((v) => v.toFixed(2));
  let value = o.value;
  const paint = () => {
    const p = (Number(input.value) - Number(input.min)) / (Number(input.max) - Number(input.min));
    input.style.setProperty('--pct', `${(p * 100).toFixed(1)}%`);
    out.textContent = f(value);
    input.setAttribute('aria-valuetext', f(value));
  };
  input.addEventListener('input', () => { value = fromPos(Number(input.value)); paint(); o.onInput?.(value); });
  const wrap = el('div', { class: 'ctl', id: `ctl-${o.id}` },
    el('label', { for: id }, el('span', {}, o.label), out), input,
    o.hint ? el('small', { id: `${id}-h` }, o.hint) : null);
  paint();
  return { el: wrap, input, get: () => value,
    set(v, fire = false) { value = v; input.value = String(toPos(v)); paint(); if (fire) o.onInput?.(v); } };
}

/** Segmented control: a radiogroup with arrow-key navigation. */
export function segmented(o) {
  let value = o.value;
  const lbl = `seg-${o.id}-${++uid}`;
  const group = el('div', { class: 'seg', role: 'radiogroup', 'aria-labelledby': lbl });
  const btns = o.options.map((op) => {
    const b = el('button', { type: 'button', role: 'radio', 'data-v': op.value }, op.label);
    b.addEventListener('click', () => api.set(op.value, true));
    b.addEventListener('keydown', (e) => {
      const d = ['ArrowRight', 'ArrowDown'].includes(e.key) ? 1 : ['ArrowLeft', 'ArrowUp'].includes(e.key) ? -1 : 0;
      if (!d) return;
      e.preventDefault();
      const i = o.options.findIndex((x) => x.value === value);
      const nx = o.options[(i + d + o.options.length) % o.options.length].value;
      api.set(nx, true);
      btns.find((x) => x.dataset.v === nx)?.focus();
    });
    group.append(b);
    return b;
  });
  const sync = () => btns.forEach((b) => { const on = b.dataset.v === value; b.setAttribute('aria-checked', String(on)); b.tabIndex = on ? 0 : -1; });
  const api = { el: el('div', { class: 'ctl', id: `ctl-${o.id}` }, el('span', { class: 'lbl', id: lbl }, o.label), group),
    get: () => value, set(v, fire = false) { value = v; sync(); if (fire) o.onChange?.(v); } };
  sync();
  return api;
}

export function toggle(o) {
  const input = el('input', { type: 'checkbox' });
  input.checked = !!o.value;
  input.addEventListener('change', () => o.onChange?.(input.checked));
  return { el: el('div', { class: 'ctl', id: `ctl-${o.id}` }, el('label', { class: 'toggle' }, input, el('span', {}, o.label))), input };
}

/** A grid of labelled numeric readouts, announced politely to screen readers. */
export function readout(keys) {
  const cells = new Map(), boxes = new Map();
  const root = el('div', { class: 'readout', role: 'status', 'aria-live': 'polite' });
  for (const { k, label } of keys) {
    const b = el('b', {}, '–');
    const box = el('div', {}, el('span', {}, label), b);
    cells.set(k, b); boxes.set(k, box); root.append(box);
  }
  return { el: root, set(k, v, tone) { const c = cells.get(k); if (c && c.textContent !== v) c.textContent = v; boxes.get(k).className = tone || ''; } };
}

export function notice(items, title = 'What to notice') {
  return el('div', { class: 'notice' }, el('h3', {}, title), el('ul', {}, items.map((t) => el('li', { html: t }))));
}

/** HiDPI canvas that follows its container's width. */
export function surface(canvas, opts, onResize) {
  const ctx = canvas.getContext('2d');
  const s = { canvas, ctx, w: 0, h: 0 };
  const fit = () => {
    const w = canvas.parentElement.clientWidth;
    if (!w) return;
    const h = Math.round(Math.min(opts.maxH ?? 9999, Math.max(opts.minH ?? 0, w * opts.aspect)));
    const dpr = Math.min(devicePixelRatio || 1, 2);
    if (w === s.w && h === s.h) return;
    s.w = w; s.h = h;
    canvas.style.height = h + 'px';
    canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    onResize?.(s);
  };
  const ro = new ResizeObserver(fit);
  ro.observe(canvas.parentElement);
  queueMicrotask(fit);
  return s;
}

/** rAF loop: runs only while playing and on-screen; starts paused under reduced motion. */
export class Loop {
  constructor(draw, watch, autoplay = true) {
    this.draw = draw; this.t = 0; this.playing = autoplay && !reducedMotion(); this.visible = true; this.dirty = true;
    this.raf = 0; this.subs = new Set();
    if (watch) new IntersectionObserver((es) => { this.visible = es.some((e) => e.isIntersecting); if (this.visible) this.kick(); }).observe(watch);
    document.addEventListener('visibilitychange', () => { if (!document.hidden) this.kick(); });
    this.kick();
  }
  kick() { if (!this.raf) { this.last = performance.now(); this.raf = requestAnimationFrame(this.tick); } }
  tick = (now) => {
    this.raf = 0;
    const dt = Math.min(0.05, (now - this.last) / 1000);
    this.last = now;
    if (!this.visible || document.hidden) return;
    if (this.playing) this.t += dt;
    if (this.playing || this.dirty) { this.dirty = false; this.draw(this.t, dt); }
    if (this.playing) this.raf = requestAnimationFrame(this.tick);
  };
  request() { this.dirty = true; this.kick(); }
  play() { this.playing = true; this.subs.forEach((f) => f(true)); this.kick(); }
  pause() { this.playing = false; this.subs.forEach((f) => f(false)); }
  toggle() { this.playing ? this.pause() : this.play(); }
}
export function playButton(loop) {
  const b = el('button', { type: 'button', class: 'btn' });
  const sync = (p) => { b.textContent = p ? '❚❚ Pause' : '▶ Play'; b.setAttribute('aria-pressed', String(p)); };
  sync(loop.playing); loop.subs.add(sync);
  b.addEventListener('click', () => loop.toggle());
  return b;
}

export function niceTicks(a, b, n = 5) {
  const span = b - a || 1, step0 = span / n, mag = 10 ** Math.floor(Math.log10(step0));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= step0) || mag * 10;
  const out = [];
  for (let v = Math.ceil(a / step) * step; v <= b + step * 1e-9; v += step) out.push(+v.toFixed(10));
  return out;
}

/** Draw axes into a box and return value→pixel mappers. */
export function axes(ctx, box, o) {
  const P = palette();
  const { x0, y0, w, h } = box;
  const X = (v) => x0 + ((v - o.x[0]) / (o.x[1] - o.x[0])) * w;
  const Y = (v) => y0 + h - ((v - o.y[0]) / (o.y[1] - o.y[0])) * h;
  ctx.save();
  ctx.font = '11px "JetBrains Mono", ui-monospace, monospace';
  ctx.strokeStyle = P.grid; ctx.fillStyle = P.muted; ctx.lineWidth = 1;
  ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  for (const t of niceTicks(o.x[0], o.x[1], o.nx ?? 6)) {
    ctx.beginPath(); ctx.moveTo(X(t), y0); ctx.lineTo(X(t), y0 + h); ctx.stroke();
    ctx.fillText((o.fx ?? String)(t), X(t), y0 + h + 4);
  }
  ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
  for (const t of niceTicks(o.y[0], o.y[1], o.ny ?? 4)) {
    ctx.beginPath(); ctx.moveTo(x0, Y(t)); ctx.lineTo(x0 + w, Y(t)); ctx.stroke();
    ctx.fillText((o.fy ?? String)(t), x0 - 6, Y(t));
  }
  ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
  if (o.xlabel) ctx.fillText(o.xlabel, x0 + w / 2, y0 + h + 32);
  if (o.ylabel) { ctx.translate(12, y0 + h / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'top'; ctx.fillText(o.ylabel, 0, 0); }
  ctx.restore();
  return { X, Y };
}

/** Deterministic PRNG + Gaussian, so demos look the same on every load. */
export function rng(seed = 1) {
  let a = seed >>> 0;
  const u = () => { a |= 0; a = (a + 0x6d2b79f5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
  const g = () => Math.sqrt(-2 * Math.log(u() || 1e-12)) * Math.cos(2 * Math.PI * u());
  return { u, g };
}

export function median(a) { const s = a.filter(Number.isFinite).sort((x, y) => x - y); return s.length ? s[s.length >> 1] : NaN; }

/** Twinkling background stars (paused under reduced motion). */
export function starfield() {
  const c = document.getElementById('stars');
  if (!c) return;
  const x = c.getContext('2d');
  let s = [];
  const r = () => {
    c.width = innerWidth * devicePixelRatio; c.height = innerHeight * devicePixelRatio;
    s = Array.from({ length: Math.round((innerWidth * innerHeight) / 5000) }, () => ({ x: Math.random() * c.width, y: Math.random() * c.height,
      r: Math.random() * 1.2 * devicePixelRatio + 0.2, a: Math.random() * 0.8 + 0.2, p: Math.random() * 6.28 }));
  };
  const d = (t) => {
    x.clearRect(0, 0, c.width, c.height);
    for (const q of s) { x.fillStyle = `rgba(220,230,255,${q.a * (0.65 + 0.35 * Math.sin(t / 1400 + q.p))})`; x.beginPath(); x.arc(q.x, q.y, q.r, 0, 6.283); x.fill(); }
    if (!reducedMotion() && !document.hidden) requestAnimationFrame(d);
  };
  addEventListener('resize', r); r(); requestAnimationFrame(d);
  document.addEventListener('visibilitychange', () => { if (!document.hidden && !reducedMotion()) requestAnimationFrame(d); });
}
