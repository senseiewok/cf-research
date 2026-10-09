// The website's layer over the CFTR model on cf/cftr/: the controls panel (a column beside the picture on a wide screen, a bottom sheet on a phone),
// touch gestures on the picture and the walk by touch (both from cftr-controls.js), and the variant box. It uses the viewer's API only; nothing is stored and nothing is sent anywhere.
import { ready } from './cftr-viewer.js';
import { parseVariant } from './cftr-variant.js';
import { sheet, touch, walkControls } from './cftr-controls.js';

const NAMES = { A: 'alanine', R: 'arginine', N: 'asparagine', D: 'aspartate', C: 'cysteine', Q: 'glutamine', E: 'glutamate', G: 'glycine', H: 'histidine', I: 'isoleucine',
                L: 'leucine', K: 'lysine', M: 'methionine', F: 'phenylalanine', P: 'proline', S: 'serine', T: 'threonine', W: 'tryptophan', Y: 'tyrosine', V: 'valine' };
const ONLY = 'This page can only place protein-level changes by residue number, such as F508del, G551D or p.Trp1282Ter, or a residue number from 1 to 1480. The residue look-up in this panel takes a number too.';
export const MESSAGES = {
  empty: 'Type one variant name or a residue number.',
  many: 'Please type one variant at a time, with no spaces.',
  range: 'Residue numbers here run from 1 to 1480.',
  dna: 'Names written at the DNA level, or inside or next to an intron (such as c.1521_1523delCTT or 621+1G>T), cannot be placed here. ' + ONLY,
  fs: 'A frameshift cannot be placed exactly here. ' + ONLY,
  unread: 'This page could not read that name. ' + ONLY,
};
const $ = (id) => document.getElementById(id);
const aa = (l) => l + ' (' + NAMES[l] + ')';

const api = await ready;
const root = $('cftr-viewer');
if (api && root) {
  const { data, length } = api, seq = data.uniprot.sequence;
  const placed = data.structures.map(st => { const m = new Map(); for (const g of st.segments) for (let i = 0; i < g.xyz.length / 3; i++) m.set(g.start + i, g.xyz.slice(3 * i, 3 * i + 3)); return m; });
  const view = (v) => api.setView({ ...v, phi: v.phi === undefined ? undefined : Math.max(0.2, Math.min(Math.PI - 0.2, v.phi)) });

  // ---- the panel (a bottom sheet below 62rem, opened by the Controls button), touch gestures on the picture and the walk by touch: cftr-controls.js
  sheet({ root, panel: $('cftr-panel'), open: $('cftr-sheet-open'), close: $('cftr-sheet-close'), title: $('cftr-panel-title'), narrow: matchMedia('(max-width: 61.99rem)') });
  touch(api);
  walkControls(api);

  // ---- the variant box: parse, check the reference letters, mark the residue and turn the picture towards it. Nothing typed reaches the address or storage.
  const form = $('cftr-variant'), input = $('cftr-variant-in'), out = $('cftr-variant-out');
  let shown = null, shownIn = null;
  const partOf = (n) => {
    const i = data.domains.filter(d => d.type === 'Domain').findIndex(d => n >= d.start && n <= d.end);
    if (i >= 0) return api.parts[i].label.toLowerCase() + ' (residues ' + api.parts[i].start + ' to ' + api.parts[i].end + ')';
    const r = data.domains.find(d => d.type === 'Region' && n >= d.start && n <= d.end);
    return r ? 'outside the four marked parts, in the stretch UniProt marks as "' + r.name + '" (' + r.start + ' to ' + r.end + ')' : 'outside the four marked parts';
  };
  const li = (t) => { const e = document.createElement('li'); e.textContent = t; return e; };
  function render() {
    if (!shown) return;
    const { start, end } = shown, idx = data.structures.findIndex(s => s.id === api.state().structure), st = data.structures[idx];
    const ns = []; for (let n = start; n <= end; n++) ns.push(n);
    const h = document.createElement('h4');
    h.textContent = (start === end ? 'Residue ' + start : 'Residues ' + start + ' to ' + end) + ' · ' + ns.map(n => aa(seq[n - 1])).join(', ');
    const ul = document.createElement('ul');
    ul.append(li('In the UniProt reference sequence: ' + ns.map(n => seq[n - 1] + ' at ' + n).join(', ') + '.'));
    ul.append(li('Part: ' + [...new Set(ns.map(partOf))].join('; ') + '.'));
    ul.append(li(st.id + ', the structure shown: ' + ns.map(n => n + (placed[idx].has(n) ? ' is placed' : ' is not placed')).join(', ') + '.'));
    if (end > start) ul.append(li('The marker is on ' + start + ', the first of them.'));
    out.replaceChildren(h, ul);
    shownIn = st.id;
  }
  function place(text) {
    const r = parseVariant(text, length);
    shown = null;
    if (!r.ok) { out.replaceChildren(Object.assign(document.createElement('p'), { textContent: MESSAGES[r.why] })); return r; }
    const wrong = r.ref.map((l, k) => [l, k ? r.end : r.start]).find(([l, n]) => l !== seq[n - 1]);
    if (wrong) {
      out.replaceChildren(Object.assign(document.createElement('p'), { textContent: 'The reference has ' + aa(seq[wrong[1] - 1]) + ' at ' + wrong[1] + ', not ' + aa(wrong[0]) + '. Nothing was placed.' }));
      return { ...r, ok: false, why: 'mismatch' };
    }
    shown = r;
    api.setWalk(r.start);
    if (api.state().playing) $('cftr-play')?.click();
    const idx = data.structures.findIndex(s => s.id === api.state().structure);
    let p = null;
    for (let k = 0; !p && k < 60; k++) p = placed[idx].get(r.start - k) || placed[idx].get(r.start + k);
    if (p) view({ theta: Math.atan2(p[1], p[0]), phi: Math.acos(p[2] / (Math.hypot(...p) || 1)) });
    render();
    return r;
  }
  form?.addEventListener('submit', (e) => { e.preventDefault(); place(input.value); });
  for (const b of document.querySelectorAll('[data-variant]')) b.addEventListener('click', () => { input.value = b.dataset.variant; place(input.value); });
  api.onChange((s) => { if (shown && s.structure !== shownIn) render(); });
}
