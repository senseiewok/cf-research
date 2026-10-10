# Trial endpoint atlas: offline tools

Build step 1 of the design in `proposals/2026-10-10-trial-endpoint-atlas.md` (board row T-0130): what interventional cystic fibrosis trials register as their primary outcome, counted from one dated, hashed snapshot of ClinicalTrials.gov.

**Nothing here has been run against the real registry.** Every test runs on synthetic fixtures with no network. The lexicon is a second draft (0.2.1-draft) for review. No model is called anywhere in this folder.

## Files

| File | What it does |
| --- | --- |
| `fetch_snapshot.py` | Fetches the snapshot through the evidence skill's `Client` (catalog gate, 1 request per second, 200-request ceiling, no redirects, honest user agent). It records `dataTimestamp`, the query strings verbatim and the request count, keeps every page, and writes `manifest.json` with the sha256 of each file and of the whole manifest. It stops (no manifest, `INCOMPLETE.txt` left) on any unexpected answer, including a received count that differs from `totalCount`, and on any unexpected error from the client. It takes the contact only from the environment variable `EVIDENCE_CONTACT` and does not run without one. Anything it prints or writes about an error has the contact and any address-shaped text replaced by `[contact]` |
| `snapshot.py` | Loads a snapshot and refuses it unless: the manifest hash matches; the manifest's file list equals the files on disk exactly, as spelled; every file matches its hash; every page a route names is listed once and lies inside the folder; route names are unique; no study appears twice in a route; and each route's page counts equal its `studies_received` and `total_count`. Every NCT id must be valid. Turns studies into plain records and outcome entries (`NCT…:P1`, …; `:P0` for a study with no primary outcome) |
| `scope.py` | Inclusion rules and flags: interventional only; condition list names CF; exclusions X1 to X6, each with a printed count and reason; "CF only" or "CF among others"; planned or actual start; the difference between the two retrieval routes |
| `lexicon.json` | The versioned rule lexicon: two-level taxonomy, one regular expression per rule with its own examples, safety subtypes, the composite rule, the not-stated rules, time-frame buckets |
| `lexicon.py` | Applies the lexicon. Each tag records the rule id and the exact matched span. Before matching, every run of whitespace in the text is collapsed to one space; spans are mapped back to the original text. A pattern with an inline `(?x)`, an escaped space or a nested quantifier is refused at load. A time frame longer than 2,000 characters is "unparseable". `python lexicon.py try "text"` shows what matches |
| `check_atlas.py` | Recomputes every count and fails on drift; rejects model tags whose quote is not an exact substring; prints the metrics on the frozen set; runs the negative controls; applies the publication stop rules |
| `make_labelling_sheet.py` | Draws the frozen set (stratified, seeded) and writes a blind sheet (wording only) and a sealed key |
| `synthetic_fixtures.py` | Synthetic studies, pages and a fake client for the tests and the default negative controls. Not registry data |
| `test_*.py` | Unit tests, one file per tool |

## Run the tests

```
cd tools/trial_atlas
python -m unittest -v
```

Standard library only. Three cases in `test_fetch_snapshot.py` use the evidence skill's real `Client`. They skip unless `cf-skills` sits next to this repository and `requests` is installed (`tools/sources/requirements.txt`). In CI they are always skipped, because the workflow does not check out `cf-skills`; that is why they were run once by hand, with the real client from the sibling checkout.

## The order of work

1. Set `EVIDENCE_CONTACT` in the environment (there is no `--contact` option, so the address never sits on a command line), then `fetch_snapshot.py --dry-run`, then the real run. **A real run is a network call and needs a person's approval and a contact address (design, build step 2; T-0117).** The default output folder is under `sources/downloads/`, which git ignores.
2. `scope.py SNAPSHOT` and `lexicon.py tag SNAPSHOT --out tags.json`.
3. `make_labelling_sheet.py --snapshot SNAPSHOT --tags tags.json --sheet SHEETDIR/sheet.csv --key KEYDIR/key.json --seed N`. The seed is required and should be chosen and written down when the draw is made. The key may not go into the sheet's folder or any folder inside it. **The person who labels must not run this command and must not see the key:** the tool reads rule results to draw the sample, and `--report` prints counts that come from them. The sheet is UTF-8 with a signature; label it in a spreadsheet and save it as "CSV UTF-8".
4. Later: model tags for the unclassified remainder, then `check_atlas.py --write-counts counts.json` and the full check with `--frozen-sheet`, `--frozen-key`, `--negative-controls`, `--canary`, `--planted`, `--control-snapshot` and `--explained`.

`check_atlas.py` exits 0 when everything passes, 1 on drift or a failed control, 2 on a usage error (including a malformed tags or model-tags file) and 3 when integrity passed but a publication stop rule tripped. `check_atlas.py --negative-controls` with no snapshot runs only the controls on the synthetic fixtures: its last line says "NOT THE GATE: controls only" and it exits 4 unless `--controls-only` is given, so it can never be mistaken for a passing gate.

A canary counts only when its entry is in the snapshot's scope, was left unclassified by the rules, and the tag is rejected for its quote. A canary whose entry is missing or already decided is an error, and the control fails.

Publication is blocked when:

- S0: the snapshot or the control snapshot carries the synthetic marker (a `_synthetic` key on a page or in `version.json`, or an `apiVersion` naming SYNTHETIC). Synthetic data can never pass the gate from a command line; the unit tests reach a passing gate only through a function argument no option sets
- S1: there is no frozen set; the frozen set does not belong to this snapshot and lexicon (the key's snapshot and lexicon hashes, and every sheet row's wording against this snapshot's text for that entry); fewer than 50 labelled rows (the design's number); any unlabelled row; a shown class with no labelled positive; or a shown class below precision 0.85 or recall 0.80 where 10 or more rows estimate it
- S2: "other" plus unclassified exceed 15% of entries. This is stricter than the design, which names "other" alone; unclassified entries are counted too until a model or a person has placed them
- S3: there is no recall route, or the routes differ by more than 10% unexplained
- S4: the controls were not run, or ran only on the synthetic fixtures (the gate needs `--canary`, `--planted` and `--control-snapshot`), or the control snapshot holds fewer than 50 outcome entries (a reviewer's number, not the design's)

The thresholds are the reviewer's judgement from the design, not a standard. They can change before the frozen set is labelled, not after.

## To confirm on the first real run

Each of these is configuration in `fetch_snapshot.DEFAULTS` or `ROUTES`, marked "to verify on the first real run":

- the largest `pageSize` (the documentation does not state it; default 100)
- that `pageToken` is the request key for the next page and `nextPageToken` the response key
- that `countTotal=true` returns `totalCount`, and whether it comes on the first page only (later pages may omit it; a different value stops the run)
- the query parameters `query.cond` and `query.term`, and the term syntax `"cystic fibrosis" OR mucoviscidosis`
- the `fields` piece names (the ones not already used by the evidence skill's provider are guesses by analogy)
- the `version` endpoint and its `dataTimestamp` key
- whether the answers fit the evidence client's 5 MB body cap at the chosen page size
- whether a data refresh during a run changes `totalCount` (the run stops if it does)

Pages are stored as the parsed JSON written back out (the client returns parsed data, not bytes), so the hashes cover what was parsed. Values and key order are kept; whitespace is not. The manifest hash is a consistency check, not a signature: anyone who can rewrite the folder can rewrite the manifest.

## Scope and the registry's terms

Scope, from the design: interventional studies whose condition list names cystic fibrosis, all statuses and dates; expanded access and observational studies out; one registry, one query, one date.

The registry's Terms and Conditions (last updated 2023-01-31, as summarised in the design) ask anyone who publishes or distributes the data to:

- name the source as ClinicalTrials.gov
- show the date the registry processed the data (`dataTimestamp`)
- describe every change made to the content
- never present the database as anything other than a United States Government database

The terms also ask that the data be kept current. A frozen snapshot cannot be, so the page must show its date prominently. Some text may belong to third parties, so quote only what each count needs. The manifest repeats these points.

## Lexicon decisions to review

- **A class was added: `exercise_capacity`** (six-minute walk, peak oxygen uptake, exercise capacity or tolerance). The design's taxonomy lacks it. The maintainer must confirm it before the frozen set is labelled.
- **Exacerbations need a CF or pulmonary context** (pulmonary, respiratory, protocol-defined, CF, PEx). A bare "exacerbations" and asthma, COPD or ABPA exacerbations stay unclassified. The non-CF negative control watches this class.
- **Liver blood tests are safety (laboratory) only:** liver function tests, liver enzymes, ALT, AST and bilirubin. The liver class is for liver-disease wording such as stiffness, steatosis, cirrhosis or CFLD.
- **DXA and DEXA are not imaging.** Body composition is tagged nutrition and growth through its own words; bone density has no class.
- **"FEV" or "forced expiratory volume" without the one-second timing** is other spirometry, never FEV1.
- **From the description** (read only when the measure gives no class), only the earliest match is kept and marked `from_description`. Generic safety wording never decides a class from there.
- **The time frame's largest written number decides the bucket, with no day-1 adjustment.** "Day 28" is up to 4 weeks; "Day 29" (often a four-week visit when day 1 is the first dose) falls in "over 4 weeks to 6 months".

## What the lexicon does not cover yet

- Negation and hedging: "FEV1 was not measured" is still tagged FEV1.
- Wording in any language but English.
- Measures with no class in the taxonomy, such as FeNO, cough counts, sleep, bone density, hearing and kidney safety tests. They stay unclassified for a model or a person.
- Ambiguous acronyms left out on purpose: CAT, FIS, NE.
- `multi_class` does not tell one instrument that spans two classes (SNOT-22: sinus and a named instrument) from two measures in one entry.
- Patient survival and graft survival are not told apart.
- Drug levels are told from nutrient and biomarker levels by a fixed list of exclusions (vitamins, glucose, proteins, cytokines and the like); a level not on the list may be read as pharmacokinetics.

## Not built in this step

- The model-assisted route itself. Only its verifier is here.
- The challenger review and its disagreement rate.
- The 20-row relabelling two weeks after the first labels.
- A command-line route for fetching the matched asthma or COPD control snapshot. `fetch_snapshot.run(route_params=...)` supports one, and the negative controls default to a synthetic control (which does not satisfy the gate).
- The site.
