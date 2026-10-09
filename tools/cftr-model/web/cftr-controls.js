// Shared by the model page (cftr-page.js) and the home page hero (cftr-home.js): a bottom sheet, touch gestures on the picture and the walk by touch.
// Each does nothing for an element a page lacks. It uses the viewer's API only; nothing is stored or sent.
const $ = (id) => document.getElementById(id);
const reduced = () => matchMedia('(prefers-reduced-motion: reduce)').matches;
const on = (t, v, f) => v.split(' ').map(x => t.addEventListener(x, f)), mk = (t, c) => Object.assign(document.createElement(t), { className: c });
const view = (api, v) => api.setView({ ...v, phi: v.phi === undefined ? undefined : Math.max(0.2, Math.min(Math.PI - 0.2, v.phi)) });

// ---- the sheet: `root` carries is-open; a closed panel is inert (out of the tab order) while `narrow` matches, or always. Opening puts `top` at the
// top of the screen while `scroll` matches, or always, the sheet at its own top and focus on the title; closing returns focus to the button if it was
// inside. No trap. A finger sets Turn or Tilt where it lands and scrolls nothing (on an iPhone it scrolled the sheet).
export function sheet({ root, panel, open, close, title, top = root, narrow = null, scroll = null }) {
  if (!root || !panel || !open) return null;
  const isOpen = () => root.classList.contains('is-open');
  const inert = () => { panel.inert = (!narrow || narrow.matches) && !isOpen(); };
  for (const r of panel.querySelectorAll('#cftr-turn,#cftr-tilt')) {
    on(r, 'touchstart', (e) => e.preventDefault());
    on(r, 'pointerdown pointermove', (e) => { const b = r.getBoundingClientRect();   // a touch stays on its first element
      if (e.pointerType == 'touch') { r.value = +r.min + (r.max - r.min) * (e.clientX - b.left - 12) / (b.width - 24); r.dispatchEvent(new Event('input', { bubbles: true })); } });
  }
  function set(o) {
    root.classList.toggle('is-open', o);
    inert();
    open.setAttribute('aria-expanded', String(o));
    if (o) {
      panel.scrollTop = 0;
      if (!scroll || scroll.matches) window.scrollTo({ top: scrollY + top.getBoundingClientRect().top, behavior: reduced() ? 'instant' : 'smooth' });   // the picture at the top, the sheet below it
      title?.focus({ preventScroll: true });
    } else if (panel.contains(document.activeElement)) open.focus({ preventScroll: true });
  }
  inert();
  narrow?.addEventListener('change', inert);
  open.addEventListener('click', () => set(!isOpen()));
  close?.addEventListener('click', () => set(false));
  document.addEventListener('keydown', (e) => { if (e.key === 'Escape' && isOpen() && !e.defaultPrevented && !e.target.closest('[role=dialog]')) set(false); });
  return set;
}

// ---- touch: one finger turns (and tilts, unless `tilt` is false), two fingers pinch to zoom, a double tap is the Reset view button.
// Mouse, wheel and keys stay with the viewer. With tilt off, the page's own touch-action (pan-y) keeps a vertical swipe for scrolling the page.
// A touch that starts on a control or the sheet is left to it.
export function touch(api, { tilt = true } = {}) {
  const canvas = api.canvas, stage = canvas.parentElement, pts = new Map();
  let base = null, lastTap = null, tap = null;
  const rebase = () => {
    const p = [...pts.values()], s = api.state();
    base = p.length >= 2 ? { d: Math.hypot(p[0].x - p[1].x, p[0].y - p[1].y) || 1, zoom: s.zoom } : null;
  };
  stage.addEventListener('pointerdown', (e) => {
    if (e.pointerType !== 'touch' || e.target.closest('#cftr-panel,[role=dialog],button,input,a,select,textarea,summary')) return;
    e.stopPropagation();
    pts.set(e.pointerId, { x: e.clientX, y: e.clientY });
    try { canvas.setPointerCapture(e.pointerId); } catch (x) { /* not capturable */ }
    tap = pts.size === 1 ? { x: e.clientX, y: e.clientY, t: e.timeStamp } : null;
    rebase();
  }, true);
  stage.addEventListener('pointermove', (e) => {
    const p = pts.get(e.pointerId);
    if (e.pointerType !== 'touch' || !p) return;
    e.stopPropagation();
    const dx = e.clientX - p.x, dy = e.clientY - p.y;
    p.x = e.clientX; p.y = e.clientY;
    if (tap && Math.hypot(e.clientX - tap.x, e.clientY - tap.y) > 10) tap = null;
    const s = api.state();
    if (pts.size === 1) view(api, { theta: s.theta + dx * 0.01, phi: tilt ? s.phi + dy * 0.01 : undefined });
    else if (base) { const q = [...pts.values()]; view(api, { zoom: Math.max(0.5, Math.min(4, base.zoom * Math.hypot(q[0].x - q[1].x, q[0].y - q[1].y) / base.d)) }); }
  }, true);
  const end = (e) => {
    if (e.pointerType !== 'touch' || !pts.has(e.pointerId)) return;
    e.stopPropagation();
    pts.delete(e.pointerId);
    try { canvas.releasePointerCapture(e.pointerId); } catch (x) { /* already released */ }
    if (e.type === 'pointerup' && tap && e.timeStamp - tap.t < 300) {
      if (lastTap && tap.t - lastTap.t < 350 && Math.hypot(tap.x - lastTap.x, tap.y - lastTap.y) < 30) { $('cftr-reset')?.click(); lastTap = null; } else lastTap = tap;
    } else lastTap = null;
    tap = null;
    rebase();
  };
  stage.addEventListener('pointerup', end, true);
  stage.addEventListener('pointercancel', end, true);
}

// ---- the walk by touch: step buttons that repeat while held, a strip to drag or tap, the number above the thumb under a finger.
// The step buttons go after `after` (default: the ruler's box). Nothing is built without the walk slider.
export function walkControls(api, { after = null } = {}) {
  const walk = $('cftr-walk'), strip = $('cftr-seq'), length = api.length;
  if (!walk) return null;
  const row = mk('div', 'cftr-nudge'), tip = mk('span', 'cftr-bubble'), at = () => api.state().walk, stop = () => clearTimeout(t);
  const go = (n) => { walk.value = Math.max(1, Math.min(length, n)); walk.dispatchEvent(new Event('input', { bubbles: true })); };
  let t, held, drag;
  for (const d of [-10, -1, 1, 10]) {
    const a = Math.abs(d), step = () => go((at() ?? 508 - d) + d), rep = (ms) => t = setTimeout(() => { held = 1; step(); rep(80); }, ms);   // no marker: 508
    const b = Object.assign(mk('button', 'cftr-btn'), { type: 'button', textContent: (d < 0 ? '−' : '+') + a, ariaLabel: (d < 0 ? 'Back ' : 'Forward ') + a + ' residue' + (a > 1 ? 's' : '') });
    on(b, 'click', (e) => { if (!e.detail || !held) step(); held = 0; });   // a tap or a key: one step; the click that ends a hold: none
    on(b, 'pointerdown', (e) => { held = 0; stop(); if (e.isPrimary && !e.button) rep(400); });
    on(b, 'pointerup pointercancel pointerleave blur', stop);
    row.append(b);
  }
  on(document, 'visibilitychange', stop);
  (after || walk.parentElement).after(row);
  walk.after(tip);
  tip.hidden = true;
  on(walk, 'pointerdown input', (e) => { if (e.pointerType == 'touch') tip.hidden = false; tip.textContent = walk.value; tip.style.left = 14 + (walk.value - 1) / (length - 1) * (walk.clientWidth - 28) + 'px'; });
  on(walk, 'pointerup pointercancel change blur', () => { tip.hidden = true; });
  if (strip) {
    on(strip, 'pointerdown', (e) => { if (e.isPrimary && !e.button) { drag = { x: e.clientX, n: at() ?? 508, w: (strip.clientWidth + 2) / strip.children.length, c: e.target.closest('span') }; strip.setPointerCapture(e.pointerId); } });
    on(strip, 'pointermove', (e) => { const k = drag && Math.round((drag.x - e.clientX) / drag.w); if (k || drag && !drag.c) { drag.c = null; go(drag.n + k); } });   // the letters follow the finger
    on(strip, 'pointerup pointercancel', (e) => { if (e.type == 'pointerup' && drag?.c) go(parseInt(drag.c.title)); drag = null; });   // a tap marks the letter tapped
  }
  return row;
}
