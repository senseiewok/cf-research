# Trial endpoint atlas: offline tools

Build step 1 of the design in `proposals/2026-10-10-trial-endpoint-atlas.md` (board row T-0130): what interventional cystic fibrosis trials register as their primary outcome, counted from one dated, hashed snapshot of ClinicalTrials.gov.

**Nothing here has been run against the real registry.** Every test runs on synthetic fixtures with no network. The lexicon is a first draft for review. No model is called anywhere in this folder.

## Files

| File | What it does |
| --- | --- |
| `fetch_snapshot.py` | Fetches the snapshot through the evidence skill's `Client` (catalog gate, 1 request per second, 200-request ceiling, no redirects, honest user agent). It records `dataTimestamp`, the query strings verbatim and the request count, keeps every page, and writes `manifest.json` with the sha256 of each file and of the whole snapshot. It stops (no manifest, `INCOMPLETE.txt` left) on any unexpected answer, including a received count that differs from `totalCount`. It does not run without a contact |
| `snapshot.py` | Loads a snapshot, re-checks every hash, and turns studies into plain records and outcome entries (`NCT…:P1`, …; `:P0` for a study with no primary outcome) |
| `scope.py` | Inclusion rules and flags: interventional only; condition list names CF; exclusions X1 to X6, each with a printed count and reason; "CF only" or "CF among others"; planned or actual start; the difference between the two retrieval routes |
| `lexicon.json` | The versioned rule lexicon: two-level taxonomy, one regular expression per rule with its own examples, safety subtypes, the composite rule, the not-stated rules, time-frame buckets |
| `lexicon.py` | Applies the lexicon. Each tag records the rule id and the exact matched span. `python lexicon.py try "text"` shows what matches |
| `check_atlas.py` | Recomputes every count and fails on drift; rejects model tags whose quote is not an exact substring; prints the metrics on the frozen set; runs the negative controls; applies the publication stop rules |
| `make_labelling_sheet.py` | Draws the frozen set (stratified, seeded) and writes a blind sheet (wording only) and a sealed key |
| `synthetic_fixtures.py` | Synthetic studies, pages and a fake client for the tests and the default negative controls. Not registry data |
| `test_*.py` | Unit tests, one file per tool |

## Run the tests

```
cd tools/trial_atlas
python -m unittest -v
```

Standard library only. Three cases in `test_fetch_snapshot.py` use the evidence skill's real `Client`; they skip unless `cf-skills` sits next to this repository and `requests` is installed (`tools/sources/requirements.txt`).

## The order of work

1. `fetch_snapshot.py --contact ADDRESS --dry-run`, then the real run. **A real run is a network call and needs a person's approval and a contact address (design, build step 2; T-0117).** The default output folder is under `sources/downloads/`, which git ignores.
2. `scope.py SNAPSHOT` and `lexicon.py tag SNAPSHOT --out tags.json`.
3. `make_labelling_sheet.py --snapshot SNAPSHOT --tags tags.json --sheet sheet.csv --key key.json`. A person labels the sheet before seeing any rule or model result.
4. Later: model tags for the unclassified remainder, then `check_atlas.py --write-counts counts.json` and the full check with `--frozen-sheet`, `--frozen-key`, `--negative-controls` and `--explained`.

`check_atlas.py` exits 0 when everything passes, 1 on drift or a failed control, 2 on a usage error and 3 when integrity passed but a publication stop rule tripped. Without a frozen set, a recall route or the controls, publication is blocked. The thresholds (precision 0.85, recall 0.80, other 15%, routes 10%, 10 positives before estimating) are the reviewer's judgement from the design, not a standard. They can change before the frozen set is labelled, not after.

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

Pages are stored as the parsed JSON written back out (the client returns parsed data, not bytes), so the hashes cover what was parsed. Values and key order are kept; whitespace is not.

## Scope and the registry's terms

Scope, from the design: interventional studies whose condition list names cystic fibrosis, all statuses and dates; expanded access and observational studies out; one registry, one query, one date.

The registry's Terms and Conditions (last updated 2023-01-31, as summarised in the design) ask anyone who publishes or distributes the data to:

- name the source as ClinicalTrials.gov
- show the date the registry processed the data (`dataTimestamp`)
- describe every change made to the content
- never present the database as anything other than a United States Government database

The terms also ask that the data be kept current. A frozen snapshot cannot be, so the page must show its date prominently. Some text may belong to third parties, so quote only what each count needs. The manifest repeats these points.

## What the lexicon does not cover yet

- Negation and hedging: "FEV1 was not measured" is still tagged FEV1.
- Wording in any language but English.
- Measures with no class in the taxonomy, such as exercise capacity (6-minute walk), FeNO, cough counts, sleep, bone density (DXA is tagged imaging only), hearing and kidney safety tests. They stay unclassified for a model or a person.
- Ambiguous acronyms left out on purpose: CAT, ALT, AST, FIS, NE.
- Exacerbations of other diseases count as exacerbations (as in the COPD control).
- `multi_class` does not tell one instrument that spans two classes (SNOT-22: sinus and a named instrument) from two measures in one entry.
- Liver function tests are tagged both liver and safety (laboratory).
- Patient survival and graft survival are not told apart.

## Not built in this step

- The model-assisted route itself. Only its verifier is here.
- The challenger review and its disagreement rate.
- The 20-row relabelling two weeks after the first labels.
- A command-line route for fetching the matched asthma or COPD control snapshot. `fetch_snapshot.run(route_params=...)` supports one, and the negative controls default to a synthetic control.
- The site.
