// The home page hero's layer over the CFTR model: a Controls button on the picture opens a bottom sheet with the walk's step buttons and the view controls,
// two fingers pinch to zoom and a double tap resets. One finger only turns: a vertical swipe on the picture still scrolls the page (touch-action: pan-y),
// so tilting is on the slider in the sheet. It uses the viewer's API only; nothing is stored and nothing is sent anywhere.
import { ready } from './cftr-viewer.js';
import { sheet, touch, walkControls } from './cftr-controls.js';

const $ = (id) => document.getElementById(id);
const api = await ready;
const root = $('cftr-viewer'), panel = $('cftr-panel');
if (api && root && panel) {
  sheet({ root, panel, open: $('cftr-sheet-open'), close: $('cftr-sheet-close'), title: $('cftr-panel-title'),
         top: root.closest('.cftr-hero-frame') || api.canvas, scroll: matchMedia('(max-width: 900px)') });   // a wide screen: the sheet opens over the text column, no scroll
  touch(api, { tilt: false });
  const row = walkControls(api, { after: panel.querySelector('.cftr-panel-head') });
  const said = $('cftr-walk-out');
  if (row && said) {               // the readout sits under the ruler, behind the open sheet: a copy beside the step buttons (the original is the one read aloud)
    const copy = Object.assign(document.createElement('p'), { className: 'meta cftr-home-walkout' });
    copy.setAttribute('aria-hidden', 'true');
    const sync = () => { copy.textContent = said.textContent; };
    new MutationObserver(sync).observe(said, { childList: true, characterData: true, subtree: true });
    sync();
    row.after(copy);
  }
  const full = $('cftr-full-link'), link = panel.querySelector('.cftr-home-full');
  if (full && link) link.addEventListener('click', () => { link.href = full.href; });   // the full page opens on the same view, as the hero's own link does
}
