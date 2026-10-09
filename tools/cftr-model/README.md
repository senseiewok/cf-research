# tools/cftr-model

An interactive 3D view of four published structures of the human CFTR protein, with research tools for reading them. It is an open, checkable example of what this lab does: take public data, say plainly what it shows and what it does not, and let a stranger verify every number.

**Not medical advice.** A structure is a still picture of one arrangement of a protein under one set of experimental conditions. Nothing here shows how any medicine acts, whether a channel is open or closed, or what anyone should expect. The research note [`landscape/cftr-structures.md`](../../landscape/cftr-structures.md) records what the files say and what we deliberately do not claim.

## What it shows

- Four human CFTR models from the RCSB Protein Data Bank (6MSM, 5UAK, 8EIQ, 8EJ1), drawn as a smooth tube through the measured alpha-carbon positions, coloured by UniProt's four marked parts, with flat bands where UniProt annotates a helix or a strand. The tube between two measured points is a drawing choice, and the page says so.
- Residue 508 as a rose ring and dot where a structure has it, and a dashed ring midway between 507 and 509 where the file says it was deleted. Stretches a model does not place are dotted straight connectors, not shapes.
- Switching structure is a cross-dissolve of two complete pictures at a fixed camera. Nothing is interpolated between structures, so nothing can look like the protein moving.
- A chain strip of all 1,480 positions, UniProt's own sequence letters around the marked residue, and a guided look in five steps.
- Research tools: look up a residue (what UniProt annotates, what each file does with it); measure an alpha-carbon distance; colour by UniProt's topology; mark UniProt's natural-variant positions (annotation only, with UniProt's own caveat); copy a link that restores the exact view; export a labelled PNG, the placed residues as CSV, or the view as JSON.
- A provenance table (where each number comes from and how to check it), a list of what the model leaves out, and a cite-this block, all generated from the data.

## The page layer: controls panel, touch and variant box

A few files sit on top of the viewer for a page that wants them. The lab's website uses them on its model page and on its home page; the standalone `index.html` here does not load them yet, because they need markup it does not have (the ids are in `cftr-page.js`, `cftr-home.js` and `cftr-controls.js`). Each script builds only what the page's markup allows and does nothing for an element that is missing, so loading them on a page without their markup is not an error (`tests/test_page_layer.py` checks this on the standalone page).

- **Shared controls** (`cftr-controls.js`, `cftr-controls.css`): the bottom sheet, touch gestures on the picture and the walk by touch, used by both pages below.
- **Controls panel** (model page, `cftr-page.js`). On a wide screen the controls sit in a column beside the picture, and both stay in view while the page scrolls. On a phone they are a sheet below the picture, opened by a Controls button and closed by its close button or Escape; a closed sheet is out of the tab order. With reduced motion asked for, the sheet does not slide.
- **Touch** (model page). One finger turns and tilts the picture, two fingers pinch to zoom, and a double tap is the Reset view button. Dragging on the picture does not scroll the page; the rest of the page scrolls as usual. Mouse, wheel and keys stay with the viewer.
- **Home page hero** (`cftr-home.js`, `cftr-home.css`). A Controls button on the picture opens a bottom sheet with the step buttons, the turn and tilt sliders, zoom, reset, the slow turn and a link to the full page; it closes with its Close button or Escape, and is out of the tab order while closed. Two fingers pinch to zoom and a double tap resets. One finger only turns: the picture fills much of a phone's first screen, so a vertical swipe on it scrolls the page (`touch-action: pan-y`), and tilting is on the slider. On a wide screen the sheet opens over the text column, not the picture. The variant box stays on the model page.
- **Variant box.** Type one protein-level name, such as F508del, G551D, p.Trp1282Ter or I507_F508del, or a residue number from 1 to 1480. The page checks the reference letters against UniProt's sequence, marks the residue, turns the picture towards it and says which part it falls in and whether the structure shown places it. It refuses DNA-level and intron names, frameshifts and more than one name at a time, and says why. The name typed never goes into the address, into storage or over the network; the address keeps only the residue number, as it does for any marked residue.

This shows where the position sits on one structure. It does not say what the change does or what it means for anyone.

CFTR2 data is not included: the page only links to CFTR2, whose terms bar republishing any of it without written permission (see the `cftr2` entry in [`sources/catalog.yaml`](../../sources/catalog.yaml)).

`cftr-variant.js` is pure (no page, no data, no network), so its cases can be tested on their own. `tests/test_page_layer.py` checks the three files are here and listed, that they make no network call and use no storage, and runs the parser's accepted and refused cases in Chromium. The panel and touch tests need the website's markup and run there.

## Run it

```bash
python tools/cftr-model/serve.py          # http://127.0.0.1:8080/ , this computer only; Ctrl+C stops it
```

It needs a browser with WebGL2. Without JavaScript or WebGL2 the page still reads: the facts are written out as text and tables.

## Rebuild the data

The data module `web/cftr-data.js` is generated; do not edit it by hand.

```bash
python tools/sources/fetch_sources.py --only rcsb-pdb-6msm rcsb-pdb-5uak rcsb-pdb-8eiq rcsb-pdb-8ej1 uniprot-p13569
python tools/cftr-model/build_traces.py       # reads sources/downloads/, writes web/cftr-data.js
python tools/cftr-model/make_page_blocks.py   # rewrites the generated tables in web/index.html
```

The fetcher follows the repo's network rules (robots.txt, the project user agent, one request at a time) and records each file's SHA-256 in `sources/manifest.json`, so a stranger can fetch the same bytes and get the same data. **Do not use `data.rcsb.org/rest`: its robots.txt disallows it.** The builder reads the metadata from the file headers instead.

What the builder does: keeps chain A, model 1, one alpha-carbon atom per residue; reads each file's declared differences from UniProt's sequence (an engineered change at 1371 in three of the four, the F508 deletion in two); refuses to continue if a placed residue disagrees with UniProt's sequence without its file saying so; fits each model onto 6MSM with Horn's quaternion method over the membrane-spanning residues both place (standard library only); and copies UniProt's sequence, parts, transmembrane and topology entries, ATP-binding sites, helix/strand/turn ranges and natural-variant entries.

## Test it

```bash
cd tools/cftr-model/tests
python -m unittest test_data test_page_blocks          # standard library only; test_data needs the fetched files
python -m unittest test_page_layer                     # the static checks need nothing; the browser cases need Playwright and are skipped without it
python -m unittest test_viewer test_research           # needs Playwright for Python and a Chromium (see below)
```

- `test_data`: the fit recovers a known rotation and never returns a mirror image; the numbering matches UniProt residue by residue; the data module is exactly what the builder makes from the fetched files.
- `test_page_blocks`: the generated parts of the page match the data, and the check fails when one is wrong.
- `test_page_layer`: the page layer's files are here and listed, make no network call and use no storage (and the detector fires on each kind of call), the variant parser accepts and refuses the same names as on the website, and the page scripts load on the standalone page (which lacks their panel and variant markup) without an error, with the shared step buttons, pinch and double tap working there.
- `test_viewer`, `test_research`: a real Chromium with a software GL, over HTTP with a strict Content-Security-Policy. They check the picture is drawn and coloured, residue 508 is marked and its absence shown, every control works with a real key or pointer, reduced motion draws one still frame and stops, a lost WebGL context is rebuilt, a phone needs no sideways scroll, the page reads without JavaScript, and every answer a research tool gives equals a value recomputed from the data (look-up against UniProt's lists, distances against the raw coordinates, CSV rows against the placed residues, variant counts against UniProt's).

Playwright is not a dependency of anything else in this repo. Install it in its own virtual environment outside the repo; the tests use a Chromium already on disk (`CHROME_PATH`, or a Playwright browser folder) and block every host except `127.0.0.1`.

The same browser tests run against another page that has the same element ids, which is how the lab's website checks its own copy:

```bash
CFTR_WEB_ROOT=/path/to/site/src CFTR_PAGE=cf/cftr/ CFTR_DATA_JS=/path/to/site/src/scripts/cftr-data.js python -m unittest test_viewer test_research
```

Tests that can never fail prove nothing, so the guards are also checked by breaking the code on purpose (a feature switched off, a number changed, a sentence removed) and confirming a test fails. On 2026-10-05, 44 such faults (34 in the viewer, 10 in the data builder) were each caught by at least one test; later that day the off-screen pause and the data-turn switch added eight more viewer faults, of which six were caught and two (removing one of the pause's three redundant guards) changed nothing a visitor or a test can see. A full re-run of all 52 on the final code gave 50 caught and 2 survived, no fault failed to apply (an unapplied fault counts as a failed proof), and the unmutated suite passed again at the end. The harness is not in this repository, so that count is a record, not something to re-run here; the review record is [`proposals/2026-10-05-cftr-model-review.md`](../../proposals/2026-10-05-cftr-model-review.md).

## Files

| Path | What it is |
| --- | --- |
| `web/index.html` | The standalone page: viewer, guided look, research tools, tables, provenance, limits |
| `web/cftr-viewer.js` | The viewer, raw WebGL2, no library and no network call; exports a small API the tools use |
| `web/cftr-research.js` | The research tools |
| `web/cftr-data.js` | Generated data: alpha-carbon traces, UniProt annotations, each file's declared differences |
| `web/cftr-model.css`, `web/cftr-tokens.css` | Styles; the second defines the colour tokens for the standalone page |
| `web/cftr-page.js` | The model page's layer: the controls panel and the variant box, with the sheet, touch and walk from `cftr-controls.js`; uses the viewer's API only |
| `web/cftr-variant.js` | A pure parser for one protein-level variant name; exports `parseVariant` |
| `web/cftr-page.css` | Styles for the model page's layer: the panel as a column or a bottom sheet, the touch hint, the variant box |
| `web/cftr-controls.js` | Shared by both pages: the bottom sheet, touch gestures on the picture (pinch, double tap, tilt optional), and step buttons and a draggable sequence strip for the walk along the chain; exports `sheet`, `touch` and `walkControls` |
| `web/cftr-controls.css` | Styles for the walk by touch: the step buttons, the strip, the wider touch thumb and its number |
| `web/cftr-home.js` | The home page hero's layer: the Controls button and its sheet, pinch and double tap, one-finger turn without tilt |
| `web/cftr-home.css` | Styles for the home page hero's Controls button and sheet |
| `build_traces.py` | Builds the data module from the fetched files |
| `make_page_blocks.py` | Writes the generated tables, the provenance and cite blocks, and a still vector picture |
| `serve.py` | A loopback preview server with a strict Content-Security-Policy |
| `tests/` | The tests above |

## Licences and attribution

- Code in this folder: MIT, like the rest of `tools/` ([`../LICENSE`](../LICENSE)). The research notes are CC BY 4.0.
- Structure data: RCSB Protein Data Bank archive files, released under CC0 1.0 by the wwPDB policy; please cite the entries and their papers, as listed on the page and in the note.
- Annotations: UniProt P13569, CC BY 4.0, credited on the page; UniProt asks that medical and genetic information be used for research and education and not as medical advice, and the page repeats that where variants are shown.
- The data module embeds values derived from those two sources with attribution; nothing is hidden in it.

## How it was made

Honesty about the process is part of the lab's rules. The controlling agent was Claude (Sonnet 5.5) under the maintainer's direction. A second model, Fable, designed the page and wrote the 3D renderer in an isolated copy of the files, with the test suite as its guard; its work was reviewed and merged by the controlling agent, who checked that nothing outside the copy changed. The lab's local worker model (Qwen3.8 27B) was tried first on writing the viewer and did not pass the verifier in six attempts (recorded in `cf-lab`). The captions were reviewed blind by a cloud model (Opus) and the local model against the checked facts; the cloud reviewer's 14 findings were each checked against a source and acted on. The process record is [`proposals/2026-10-05-cftr-model-review.md`](../../proposals/2026-10-05-cftr-model-review.md).

## Known limits

Draft: a person with cystic fibrosis or a carer has not yet read the wording, and the science has not had an outside expert's check. Alpha-carbon traces only; chain A of model 1 only; the four entries come from different samples and conditions; the compounds in an entry are not drawn; and CFTR2 and registry data were not read, so nothing is said about frequency or severity.
