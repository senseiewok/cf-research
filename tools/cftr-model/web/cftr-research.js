// Tools for reading the four CFTR structures, on the research view (cf/cftr/). They use the viewer's API and the same data module; nothing is sent anywhere.
//  - look up a residue: what UniProt annotates there, and what each of the four structures does with it
//  - colour the chain by UniProt's topology, or mark UniProt's annotated natural-variant positions (both are annotation, labelled as such)
//  - measure the alpha-carbon to alpha-carbon distance between two placed residues
//  - copy a link to the current view; export the placed residues (CSV), the view (JSON) or a picture with its label baked in (PNG)
// Every sentence a tool writes is built from the data; the fixed notes beside the tools are in the page's HTML.
import { ready } from './cftr-viewer.js';

const NAMES = { A: 'alanine', R: 'arginine', N: 'asparagine', D: 'aspartate', C: 'cysteine', Q: 'glutamine', E: 'glutamate', G: 'glycine', H: 'histidine', I: 'isoleucine',
                L: 'leucine', K: 'lysine', M: 'methionine', F: 'phenylalanine', P: 'proline', S: 'serine', T: 'threonine', W: 'tryptophan', Y: 'tyrosine', V: 'valine' };
const ONE = { ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G', HIS: 'H', ILE: 'I', LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S', THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V' };
const TOPOLOGY_COLOURS = { transmembrane: [0.45, 0.55, 1.0], cytoplasmic: [0.25, 0.9, 0.55], extracellular: [0.95, 0.95, 1.0], none: [0.55, 0.6, 0.66] };
const VARIANT_COLOUR = [0.85, 0.75, 1.0];

const $ = (id) => document.getElementById(id);
const make = (tag, props = {}, ...kids) => {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(props)) { if (k === 'class') e.className = v; else e.setAttribute(k, v); }
  e.append(...kids);
  return e;
};
const range = (a, b) => (a === b ? String(a) : a + ' to ' + b);

const api = await ready;
if (api) {
  const { data, length } = api;
  const u = data.uniprot;
  const placed = data.structures.map(st => {
    const at = new Map();
    for (const seg of st.segments) for (let i = 0; i < seg.xyz.length / 3; i++) at.set(seg.start + i, seg.xyz.slice(3 * i, 3 * i + 3));
    return at;
  });
  const allIds = data.structures.map(s => s.id);
  const currentIndex = () => allIds.indexOf(api.state().structure);

  // ---- what UniProt annotates at a residue (all from the data module) ----
  const within = (list, n) => list.filter(x => n >= x.start && n <= x.end);
  const variantsAt = (n) => u.variants.filter(v => n >= v.start && n <= v.end);
  const partAt = (n) => {
    const out = [];
    const domains = data.domains.filter(x => x.type === 'Domain');
    domains.forEach((x, i) => { if (n >= x.start && n <= x.end) out.push(api.parts[i].label.toLowerCase() + ' (UniProt: "' + x.name + '"), residues ' + range(x.start, x.end)); });
    for (const x of data.domains.filter(y => y.type === 'Region')) if (n >= x.start && n <= x.end) out.push('UniProt region "' + x.name + '", residues ' + range(x.start, x.end));
    return out;
  };
  const topologyAt = (n) => {
    const tm = within(u.transmembrane, n)[0];
    if (tm) return { kind: 'transmembrane', text: 'transmembrane segment, ' + tm.name.replace('Helical; Name=', 'helix ') + ' (residues ' + range(tm.start, tm.end) + ')' };
    const t = within(u.topology, n)[0];
    if (t) return { kind: t.name.toLowerCase(), text: t.name.toLowerCase() + ' side (residues ' + range(t.start, t.end) + ')' };
    return { kind: 'none', text: 'no topology entry covers this residue' };
  };
  const firstClause = (v) => (v.description.split(';')[0] || v.description).trim();
  const variantLabel = (v) => {
    const sub = v.from && v.to.length ? NAMES[v.from] + ' to ' + v.to.map(x => NAMES[x] || x).join(' or ') : 'UniProt gives no replacement residue for this entry';
    return v.id + (v.start === v.end ? ' at ' + v.start : ' spanning residues ' + range(v.start, v.end)) + ': ' + sub;
  };
  const status = (idx, n) => {
    const st = data.structures[idx];
    if (placed[idx].has(n)) return 'placed in this model';
    const d = st.differences.find(x => x.position === n && !x.file_residue);
    if (d) return 'not in the file: the file records it as "' + d.details + '"';
    let a = n, b = n;
    while (a > 1 && !placed[idx].has(a)) a--;
    while (b < length && !placed[idx].has(b)) b++;
    return 'not placed (the nearest placed residues are ' + (placed[idx].has(a) ? a : 'none before') + ' and ' + (placed[idx].has(b) ? b : 'none after') + ')';
  };

  function dossier(n) {
    const out = $('cftr-dossier');
    if (!out) return;
    const letter = u.sequence[n - 1];
    const kids = [make('h4', {}, 'Residue ' + n + ' · ' + NAMES[letter] + ' (' + letter + ')')];
    const ann = make('ul', {});
    const parts = partAt(n);
    ann.append(make('li', {}, 'Part: ' + (parts.length ? parts.join('; ') : 'outside the four marked parts')));
    ann.append(make('li', {}, 'Topology: ' + topologyAt(n).text));
    const sites = within(u.atp_sites, n);
    ann.append(make('li', {}, 'ATP-binding site annotation: ' + (sites.length ? 'yes (' + sites.map(s => range(s.start, s.end)).join(', ') + ')' : 'none at this residue')));
    const vs = variantsAt(n);
    const vli = make('li', {}, 'Annotated natural variants here: ' + vs.length);
    if (vs.length) {
      const list = make('ul', {});
      for (const v of vs) {
        const item = make('li', {}, variantLabel(v) + '. UniProt says: "' + firstClause(v) + '"');
        if (v.description.trim() !== firstClause(v)) {
          item.append(make('details', {}, make('summary', {}, 'UniProt’s full wording'), make('p', {}, v.description)));
        }
        list.append(item);
      }
      vli.append(list);
    }
    ann.append(vli);
    kids.push(make('p', { class: 'cftr-dossier-group' }, 'UniProt annotation (P13569, CC BY 4.0). This is annotation, not a clinical statement.'), ann);
    const rows = make('ul', {});
    data.structures.forEach((st, i) => {
      const li = make('li', {}, st.id + ': ' + status(i, n));
      const eng = st.differences.find(x => x.position === n && x.file_residue);
      if (eng) li.append(' — the file records ' + NAMES[ONE[eng.file_residue]] + ' here where UniProt has ' + NAMES[u.sequence[n - 1]] + ' ("' + eng.details + '").');
      rows.append(li);
    });
    kids.push(make('p', { class: 'cftr-dossier-group' }, 'From the four structures (read from their files)'), rows);
    kids.push(make('p', {}, make('button', { type: 'button', class: 'cftr-btn', id: 'cftr-dossier-mark' }, 'Mark residue ' + n + ' on the model')));
    out.replaceChildren(...kids);
    $('cftr-dossier-mark').addEventListener('click', () => api.setWalk(n));
  }
  const lookup = $('cftr-lookup');
  if (lookup) {
    const input = $('cftr-lookup-n');
    lookup.addEventListener('submit', (e) => {
      e.preventDefault();
      const n = Number(input.value);
      if (!Number.isInteger(n) || n < 1 || n > length) { $('cftr-dossier').replaceChildren(make('p', {}, 'Please enter a whole number from 1 to ' + length + '.')); return; }
      dossier(n);
    });
    input.value = '508';
  }

  // ---- distance between two placed residues, alpha carbon to alpha carbon ----
  const measureForm = $('cftr-measure');
  let measured = null;
  const measureLayer = (st, idx) => {
    if (!measured || measured.idx !== idx) return [];
    const [a, b] = measured.pts, d = measured.d, k = Math.max(3, Math.round(d / 2.5));
    const dots = [];
    for (let j = 1; j < k; j++) dots.push({ pos: [a[0] + (b[0] - a[0]) * j / k, a[1] + (b[1] - a[1]) * j / k, a[2] + (b[2] - a[2]) * j / k], size: 5, colour: [1, 1, 1], mode: 0, alpha: 0.85 });
    dots.push({ pos: a, size: 12, colour: [1, 1, 1], mode: 0 }, { pos: b, size: 12, colour: [1, 1, 1], mode: 0 });
    return dots;
  };
  if (measureForm) {
    measureForm.addEventListener('submit', (e) => {
      e.preventDefault();
      const out = $('cftr-measure-out'), from = Number($('cftr-measure-a').value), to = Number($('cftr-measure-b').value), idx = currentIndex();
      const ok = (n) => Number.isInteger(n) && n >= 1 && n <= length;
      measured = null;
      api.setLayer('measure', false);
      if (!ok(from) || !ok(to)) { out.textContent = 'Please enter two whole numbers from 1 to ' + length + '.'; return; }
      const missing = [from, to].filter(n => !placed[idx].has(n));
      if (missing.length) { out.textContent = 'Not measurable in ' + allIds[idx] + ': residue ' + missing.join(' and ') + ' ' + (missing.length > 1 ? 'are' : 'is') + ' not placed in this model.'; return; }
      const a = placed[idx].get(from), b = placed[idx].get(to), d = Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
      measured = { idx, pts: [a, b], d };
      api.setLayer('measure', true, measureLayer);
      out.textContent = allIds[idx] + ': residue ' + from + ' to residue ' + to + ' is ' + d.toFixed(1) + ' Å, alpha carbon to alpha carbon, in this model. Positions are stored to 0.1 Å. This is not a contact or a bond.';
    });
    api.onChange(() => { /* keep the line only for the structure it was measured in */ });
  }

  // ---- layers: UniProt topology colouring, UniProt variant positions ----
  const topologyBox = $('cftr-layer-topology');
  const topologyColour = (n) => TOPOLOGY_COLOURS[topologyAt(n).kind] || TOPOLOGY_COLOURS.none;
  const topologyCache = new Array(length + 2);
  const topologyFast = (n) => topologyCache[n] || (topologyCache[n] = topologyColour(n));
  function setTopology(on) {
    api.setColour(on ? topologyFast : null);
    api.setLayer('topology', on);
    const legend = $('cftr-topology-legend');
    if (legend) legend.hidden = !on;
    if (topologyBox) topologyBox.checked = on;
  }
  on_(topologyBox, 'change', () => setTopology(topologyBox.checked));

  const variantsBox = $('cftr-layer-variants');
  const singleVariants = u.variants.filter(v => v.start === v.end);
  const spanningVariants = u.variants.filter(v => v.start !== v.end);
  const variantPositions = [...new Set(singleVariants.map(v => v.start))].sort((a, b) => a - b);
  const variantMarkers = data.structures.map((st, i) => variantPositions.filter(n => placed[i].has(n)).map(n => ({ pos: placed[i].get(n), size: 9, colour: VARIANT_COLOUR, mode: 0, alpha: 0.9 })));
  const variantLayer = (st, idx) => variantMarkers[idx];
  function variantSummary() {
    const box = $('cftr-variant-summary');
    if (!box) return;
    const idx = currentIndex(), here = variantPositions.filter(n => placed[idx].has(n)).length;
    box.textContent = u.variants.length + ' UniProt entries: ' + singleVariants.length + ' at a single residue (' + variantPositions.length + ' different positions) and ' + spanningVariants.length
      + ' spanning several residues, which are not drawn (look the residues up to see them). ' + allIds[idx] + ' places ' + here + ' of the ' + variantPositions.length + ' single-residue positions; the other '
      + (variantPositions.length - here) + ' are not placed in it, so nothing is drawn for them.';
  }
  function setVariants(on) {
    api.setLayer('variants', on, on ? variantLayer : undefined);
    if (variantsBox) variantsBox.checked = on;
    variantSummary();
  }
  on_(variantsBox, 'change', () => setVariants(variantsBox.checked));
  api.onChange(() => variantSummary());
  for (const name of api.layerNames()) { if (name === 'topology') setTopology(true); if (name === 'variants') setVariants(true); }
  variantSummary();

  // ---- sharing and exports ----
  const say = (t) => { const s = $('cftr-export-status'); if (s) s.textContent = t; };
  function download(name, type, content) {
    const url = URL.createObjectURL(new Blob([content], { type }));
    const a = make('a', { href: url, download: name });
    document.body.append(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 2000);
  }
  on_($('cftr-copy-link'), 'click', async () => {
    const link = location.href;
    try { await navigator.clipboard.writeText(link); say('The address of this view was copied.'); }
    catch (e) { say('Copying was not allowed here. The address bar holds the link to this view.'); }
  });
  on_($('cftr-export-csv'), 'click', () => {
    const idx = currentIndex(), st = data.structures[idx], rows = ['structure,residue,amino_acid,x_angstrom,y_angstrom,z_angstrom,uniprot_part,uniprot_topology'];
    for (const [n, p] of [...placed[idx].entries()].sort((a, b) => a[0] - b[0])) {
      const part = data.domains.filter(d => d.type === 'Domain' && n >= d.start && n <= d.end).map(d => d.name)[0] || '';
      rows.push([st.id, n, u.sequence[n - 1], p[0], p[1], p[2], '"' + part + '"', topologyAt(n).kind].join(','));
    }
    download(st.id + '-alpha-carbons.csv', 'text/csv', rows.join('\n') + '\n');
    say('Saved ' + (rows.length - 1) + ' placed residues of ' + st.id + '.');
  });
  on_($('cftr-export-json'), 'click', () => {
    const s = api.state();
    download('cftr-view.json', 'application/json', JSON.stringify({ address: location.href, structure: s.structure, overlay: s.overlay, marked_residue: s.walk, view: { theta_radians: s.theta, phi_radians: s.phi, zoom: s.zoom }, layers: s.layers,
      data: 'RCSB PDB (CC0) and UniProt P13569 (CC BY 4.0); see the page for the files' }, null, 2) + '\n');
    say('Saved the current view.');
  });
  on_($('cftr-export-png'), 'click', () => {
    api.redrawNow();
    api.canvas.toBlob((blob) => {
      if (!blob) { say('The picture could not be saved here.'); return; }
      createImageBitmap(blob).then((bmp) => {
        const bar = 64, c = make('canvas', {}); c.width = bmp.width; c.height = bmp.height + bar;
        const g = c.getContext('2d');
        g.fillStyle = '#050508'; g.fillRect(0, 0, c.width, c.height); g.drawImage(bmp, 0, 0);
        const s = api.state();
        g.fillStyle = '#e6e2f0'; g.font = Math.max(12, Math.round(c.width / 60)) + 'px monospace';
        g.fillText('CFTR structure ' + s.structure + ' · alpha-carbon trace · a model fitted to electron-microscope data, not a photograph', 12, bmp.height + 24);
        g.fillText('Data: RCSB PDB (CC0), UniProt P13569 (CC BY 4.0) · ' + location.origin + location.pathname + '#' + api.encodeHash(), 12, bmp.height + 48);
        c.toBlob((out) => { const url = URL.createObjectURL(out); const a = make('a', { href: url, download: 'cftr-' + s.structure + '.png' }); document.body.append(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(url), 2000); say('Saved a picture with its label.'); });
      });
    }, 'image/png');
  });
  for (const n of document.querySelectorAll('[data-needs-gl]')) n.hidden = false;
}

function on_(node, type, fn) { if (node) node.addEventListener(type, fn); }
