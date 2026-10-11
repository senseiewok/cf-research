# Trial endpoint atlas: offline tools

Build step 1 of the design in `proposals/2026-10-10-trial-endpoint-atlas.md` (board row T-0130): what interventional cystic fibrosis trials register as their primary outcome, counted from one dated, hashed snapshot of ClinicalTrials.gov.

**Nothing here has been run against the real registry.** Every test runs on synthetic fixtures with no network. The lexicon is a draft (0.2.2-draft) for review. No model is called anywhere in this folder.

Words used here: a **retrieval route** is one of the searches a snapshot is fetched with (the condition search, which is the main retrieval route, and the term search, which is the recall retrieval route). The **route name** is that search's name as the manifest records it (`condition`, `term`).

## Files

| File | What it does |
| --- | --- |
| `fetch_snapshot.py` | Fetches the snapshot through the evidence skill's `Client` (catalog gate, 1 request per second, 200-request ceiling, no redirects, honest user agent). It records `dataTimestamp`, the query strings verbatim and the request count, keeps every page, and writes `manifest.json` with the sha256 of each file and of the whole manifest. It stops (no manifest, `INCOMPLETE.txt` left) on any unexpected answer, including a received count that differs from `totalCount`, and on any unexpected error from the client. It takes the contact only from the environment variable `EVIDENCE_CONTACT` and does not run without one. Anything it prints or writes about an error has the contact and any address-shaped text replaced by `[contact]` |
| `snapshot.py` | Loads a snapshot and refuses it unless: the manifest hash matches; the manifest's file list equals the files on disk exactly, as spelled; every file matches its hash; every page a retrieval route names is listed once and lies inside the folder; route names are unique; no study appears twice in a retrieval route; and each retrieval route's page counts equal its `studies_received` and `total_count`. Every NCT id must be valid. Turns studies into plain records and outcome entries (`NCT…:P1`, …; `:P0` for a study with no primary outcome) |
| `scope.py` | Inclusion rules and flags: interventional only; condition list names CF; exclusions X1 to X6, each with a printed count and reason; "CF only" or "CF among others"; planned or actual start; the difference between the two retrieval routes |
| `lexicon.json` | The versioned rule lexicon: two-level taxonomy (each class with a one-line gloss for the person labelling), one regular expression per rule with its own examples, safety subtypes, the composite rule, the not-stated rules, time-frame buckets |
| `lexicon.py` | Applies the lexicon. Each tag records the rule id and the exact matched span. Before matching, every run of whitespace in the text is collapsed to one space; spans are mapped back to the original text. A pattern with an inline `(?x)`, an escaped space or a nested quantifier is refused at load. A time frame longer than 2,000 characters is "unparseable". `lexicon.py try "text"` shows what matches; `lexicon.py classes` lists every class with its domain and gloss; `lexicon.py tag` never overwrites its output |
| `check_atlas.py` | Recomputes every count and fails on drift; rejects model tags whose quote is not an exact substring; prints the metrics on the frozen set; runs the negative controls; applies the publication stop rules; `--confirm-real-run` records that a person checked a real snapshot's configuration |
| `make_labelling_sheet.py` | Draws the frozen set (stratified, seeded) and writes a blind sheet (wording only), the registry's terms beside it, and a sealed key |
| `synthetic_fixtures.py` | Synthetic studies, pages and a fake client for the tests and the default negative controls. Not registry data |
| `test_*.py` | Unit tests, one file per tool |

## Run the tests

```
cd tools/trial_atlas
python -m unittest -v
```

Standard library only. Three cases in `test_fetch_snapshot.py` use the evidence skill's real `Client`. They skip unless `cf-skills` sits next to this repository and `requests` is installed. **They are always skipped in CI**, because the workflow does not check out `cf-skills`. They were run by hand during the third review round, from `tools/trial_atlas`. `requests` is not installed on the build machine, so a stand-in `requests` module that opens no socket was put on the path; the evidence skill's `Client` code itself was the real one, from the `cf-skills` checkout:

```
PYTHONPATH=<folder holding the stand-in requests module> python -m unittest -v test_fetch_snapshot.RealClientTest
```

The result seen was `Ran 3 tests in 0.059s` and `OK`. That is all that is claimed for them.

## The order of work

**Prerequisites:**
- Python 3.10 or newer.
- `cf-skills` cloned beside this repository, so that `../cf-skills/.claude/skills/cf-evidence-loop/scripts` exists. `EVIDENCE_SKILL_DIR` can point elsewhere.
- `pip install -r tools/sources/requirements.txt`, for `requests` and PyYAML.
- `EVIDENCE_CONTACT` set in the shell session that runs the fetch, for example `$env:EVIDENCE_CONTACT = '<contact address>'` in PowerShell. There is no `--contact` option, so the address never sits on a command line.

Every command below is run from the repository root. `SNAP` stands for the snapshot folder that step 1 printed.

**Who runs what.** The person who labels the frozen set must not run steps 3 and 4 (tagging and the draw), must not read their output, and must not see the key. Steps 3 and 4 are run by someone else, an agent or a second person, who hands over only `sheet.csv` and its `sheet.TERMS.txt`. A maintainer working alone must not run them; ask an agent to.

1. **Fetch the CF snapshot** (a network call; it needs a person's approval and a contact address: design, build step 2; T-0117). Check the plan first, then fetch:

   ```
   python tools/trial_atlas/fetch_snapshot.py --dry-run --out sources/downloads/trial_atlas/cf-snapshot-1
   python tools/trial_atlas/fetch_snapshot.py --out sources/downloads/trial_atlas/cf-snapshot-1
   ```

   - The dry run prints the planned requests, the absolute output folder and the page ceiling. The ceiling is `1 + retrieval routes x --max-pages`, and it must fit the client's 200-request budget, or the run is refused before any request.
   - The real run prints the absolute output folder at the end.
   - `--fields A,B,...` replaces the field piece names if the first real run shows one is wrong.
   - `sources/downloads/` is ignored by git.

2. **Fetch the matched non-CF control snapshot and confirm both configurations.** Name the retrieval route `condition`, because the tagger control reads the main retrieval route:

   ```
   python tools/trial_atlas/fetch_snapshot.py --route-query "condition=query.cond=asthma OR COPD" --out sources/downloads/trial_atlas/control-1
   python tools/trial_atlas/check_atlas.py --snapshot SNAP --confirm-real-run
   python tools/trial_atlas/check_atlas.py --snapshot sources/downloads/trial_atlas/control-1 --confirm-real-run
   ```

   - `--route-query NAME=KEY=VALUE` (repeatable) replaces the built-in retrieval routes. `KEY` must be a `query.*` or `filter.*` parameter.
   - Confirm only after a person has compared the fetched pages with the registry's documentation (see "To confirm on the first real run"). Confirming changes the snapshot hash, so do it before step 3.

3. **Scope and tags** (an agent or a second person, not the labeller):

   ```
   python tools/trial_atlas/scope.py SNAP
   python tools/trial_atlas/lexicon.py tag SNAP --out sources/downloads/trial_atlas/tags-1.json
   ```

4. **Draw the frozen set** (the same agent or second person). The seed is required: choose it when the draw is made and write it down. The key may not go into the sheet's folder or any folder inside it.

   ```
   python tools/trial_atlas/make_labelling_sheet.py --snapshot SNAP --tags sources/downloads/trial_atlas/tags-1.json --sheet handover/sheet.csv --key sealed/key.json --seed 4711
   ```

   - The tool also writes `handover/sheet.TERMS.txt`, the registry's terms for the text in the sheet.
   - A warning is printed when fewer rows are drawn than `--n` asks for.
   - `--report` prints how many rows each oversample took. That number comes from the rules, so the labeller must not see it.

5. **Label** (the person). Open `sheet.csv` in a spreadsheet, fill `label_classes` with class ids separated by `;`, and save it as "CSV UTF-8". `python tools/trial_atlas/lexicon.py classes` lists the class ids with a one-line gloss each.

6. **Model tags** for the unclassified remainder: not built yet. A model may propose a class only with an exact quote; the checker verifies every quote.

7. **Counts:**

   ```
   python tools/trial_atlas/check_atlas.py --snapshot SNAP --tags sources/downloads/trial_atlas/tags-1.json --model-tags model-tags.json --write-counts counts.json
   ```

8. **The gate:**

   ```
   python tools/trial_atlas/check_atlas.py --snapshot SNAP --tags sources/downloads/trial_atlas/tags-1.json --counts counts.json --model-tags model-tags.json --frozen-sheet handover/sheet.csv --frozen-key sealed/key.json --explained explained.json --negative-controls --canary canary.json --planted planted.json --control-snapshot sources/downloads/trial_atlas/control-1
   ```

### Files the gate reads

| File | Format |
| --- | --- |
| `--canary` | A JSON list of `{"entry_id": "NCT…:P1", "class": "<class id>", "quote": "<text>"}`. Each canary's entry must be in scope in SNAP and left unclassified by the rules, and its quote must NOT be an exact piece of that entry's wording. Any other kind of canary is an error |
| `--planted` | A JSON list of NCT ids of studies known not to be CF studies. They are looked for in both retrieval routes, and scope must exclude every one. An id in neither retrieval route fails the control |
| `--explained` and `--exclude` | A JSON object mapping each NCT id to a non-empty reason |
| `--model-tags` | A JSON list of `{entry_id, class, quote}` objects, as for canaries but with true quotes |
| `--control-snapshot` | A snapshot folder fetched with `--route-query condition=...` (step 2) |

A missing or unreadable input file is a usage error: it exits 2 with an `ERROR:` line.

### Exit codes and stop rules

`check_atlas.py` exits with:
- 0 when everything passes
- 1 on drift or a failed control
- 2 on a usage error
- 3 when integrity passed but a publication stop rule tripped
- 4 for `--negative-controls` with no snapshot: that run only exercises the controls on the synthetic fixtures, prints "NOT THE GATE: controls only", and exits 0 only with `--controls-only`.

Publication is blocked when:

- **S0:** the snapshot or the control snapshot carries the synthetic marker (a `_synthetic` key on a page or in `version.json`, or an `apiVersion` naming SYNTHETIC). Synthetic data can never pass the gate from a command line; the unit tests reach a passing gate only through a function argument that no option sets.
- **S1:** any of these:
  - there is no frozen set
  - the frozen set does not belong to this snapshot and lexicon. The key's snapshot and lexicon hashes are checked, and every sheet row's wording is compared with this snapshot's text for its entry
  - fewer than 50 rows are labelled (the design's number)
  - any row is unlabelled
  - a shown class has no labelled positive
  - a shown class falls below precision 0.85 or recall 0.80 where 10 or more rows estimate it.
- **S2:** "other" plus unclassified exceed 15% of entries.
- **S3:** there is no recall retrieval route, or the two retrieval routes differ by more than 10% unexplained.
- **S4:** any of these:
  - the controls were not run
  - they ran only on the synthetic fixtures (the gate needs `--canary`, `--planted` and `--control-snapshot`)
  - the control snapshot holds fewer than 50 outcome entries (a reviewer's number, not the design's).

**The gate is stricter than the design in two places:**
- S1 blocks any shown class with no labelled positive, where the design would only say "few checked".
- S2 counts "other" plus unclassified, where the design names "other" alone.

The proposal is being updated to match. The thresholds are the reviewer's judgement from the design, not a standard. They can change before the frozen set is labelled, not after.

## When a fetch stops

A fetch that stops part-way:
- writes `INCOMPLETE.txt` in its output folder, holding the reason and the number of requests made
- writes no `manifest.json`, so no other tool will read the folder
- leaves the folder in place, and prints its absolute path on stderr.

A stopped folder is never reused:
- the next run without `--out` makes a new folder named with a new UTC time stamp
- with `--out`, a different folder must be named, because no folder is ever overwritten.

Look first at:
1. the reason line of `INCOMPLETE.txt`
2. the request count on stderr
3. the last `page-NNNN.json` in the folder.

Do not retry in a loop: the client's per-process budget and cooldown rules exist to prevent that.

## To confirm on the first real run

Each of these is configuration in `fetch_snapshot.DEFAULTS` or `ROUTES`, marked "to verify on the first real run". After a person has checked them against the pages of a real run, `check_atlas.py --snapshot SNAP --confirm-real-run` sets `config_verified` in that snapshot's manifest. It is false in every manifest until then, and a synthetic snapshot cannot be confirmed.

- the largest `pageSize` (the documentation does not state it; default 100)
- that `pageToken` is the request key for the next page and `nextPageToken` the response key
- that `countTotal=true` returns `totalCount`, and whether it comes on the first page only (later pages may omit it; a different value stops the run)
- the query parameters `query.cond` and `query.term`, and the term syntax `"cystic fibrosis" OR mucoviscidosis`
- the `fields` piece names (the ones not already used by the evidence skill's provider are guesses by analogy; `--fields` overrides them)
- the `version` endpoint and its `dataTimestamp` key
- whether the answers fit the evidence client's 5 MB body cap at the chosen page size
- whether a data refresh during a run changes `totalCount` (the run stops if it does)

Pages are stored as the parsed JSON written back out (the client returns parsed data, not bytes), so the hashes cover what was parsed. Values and key order are kept; whitespace is not. The manifest hash is a consistency check, not a signature: anyone who can rewrite the folder can rewrite the manifest.

### Registry-side risks

- **A 403 is possible.** The registry may refuse an agent that is not a browser. The evidence client treats a 403 as a refusal, puts the host in cooldown, and does not retry. The run then stops with `INCOMPLETE.txt`. Do not work around it with a browser user agent: the client has no way to set one, on purpose.
- **robots.txt and the API terms disagree.** The registry's `robots.txt` disallows `/api/` for crawlers, while the catalog entry (`sources/catalog.yaml`, `clinicaltrials-gov`) records that the published data API's terms govern programmatic use, and the evidence client applies robots rules to `fetch` sources only. A person should re-read both before the first real run.
- **Unconfirmed configuration.** `config_verified` stays false until a person runs `--confirm-real-run`.

## Scope and the registry's terms

Scope, from the design: interventional studies whose condition list names cystic fibrosis, all statuses and dates; expanded access and observational studies out; one registry, one query, one date.

The registry's Terms and Conditions (last updated 2023-01-31, as summarised in the design) ask anyone who publishes or distributes the data to:

- name the source as ClinicalTrials.gov
- show the date the registry processed the data (`dataTimestamp`)
- describe every change made to the content
- never present the database as anything other than a United States Government database

The terms also ask that the data be kept current. A frozen snapshot cannot be, so the page must show its date prominently. Some text may belong to third parties, so quote only what each count needs.

**The terms travel with the data.** `counts.json` carries a `registry_terms` block, and `check_atlas.py` prints the same block at the top of every check. The drift check covers it like any other count. The block holds:
- the source and the registry's processing date
- the snapshot's fetch date
- the terms page (`https://clinicaltrials.gov/about-site/terms-conditions`, the address given at review; it was not opened during this offline build) and the last-updated date the manifest records
- every modification the lab made: the classification by the versioned lexicon, the scope exclusions with their counts, no clipping, pages re-saved after parsing
- the licence line: "The lab's licence covers its own tags, code and counts only; registry text and fields remain ClinicalTrials.gov data under its terms."
- the note that the terms apply for as long as the data are kept.

The only CSV the tools write is the labelling sheet. A CSV has no room for a header comment, so the same block is written beside it as `<sheet name>.TERMS.txt`.

## Lexicon decisions to review

- **A class was added: `exercise_capacity`** (six-minute walk, peak oxygen uptake, exercise capacity or tolerance). The proposal now lists it, marked as pending the maintainer's confirmation, which must come before the frozen set is labelled.
- **Exacerbations need a CF or pulmonary context** (pulmonary, respiratory, protocol-defined, CF, PEx). A bare "exacerbations" and asthma, COPD or ABPA exacerbations stay unclassified. The non-CF negative control watches this class.
- **Liver blood tests are safety (laboratory) only:** liver function tests, liver enzymes, ALT, AST and bilirubin. The liver class is for liver-disease wording such as stiffness, steatosis, cirrhosis or CFLD.
- **DXA and DEXA are not imaging.** Body composition is tagged nutrition and growth through its own words; bone density has no class.
- **"FEV" or "forced expiratory volume" without the one-second timing** is other spirometry, never FEV1.
- **From the description** (read only when the measure gives no class), only the earliest match is kept and marked `from_description`. Generic safety wording never decides a class from there.
- **The time frame's largest written number decides the bucket, with no day-1 adjustment.** "Day 28" is up to 4 weeks; "Day 29" (often a four-week visit when day 1 is the first dose) falls in "over 4 weeks to 6 months".
- **The glosses** (`lexicon.py classes`) are short, neutral aids for the person labelling, not text for the public site.

## What the lexicon does not cover yet

- Negation and hedging: "FEV1 was not measured" is still tagged FEV1.
- Wording in any language but English.
- Measures with no class in the taxonomy, such as FeNO, cough counts, sleep, bone density, hearing and kidney safety tests. They stay unclassified for a model or a person.
- Ambiguous acronyms left out on purpose: CAT, FIS, NE.
- `multi_class` does not tell one instrument that spans two classes (SNOT-22: sinus and a named instrument) from two measures in one entry.
- Patient survival and graft survival are not told apart.
- Drug levels are told from nutrient and biomarker levels by a fixed list of exclusions (vitamins, glucose, proteins, cytokines and the like); a level not on the list may be read as pharmacokinetics.

## Not built in this step

- The model-assisted tagging step itself. Only its quote verifier is here.
- The challenger review and its disagreement rate.
- The relabelling check: 20 rows labelled again two weeks after the first labels.
- The hand check of 20 studies from the difference between the two retrieval routes (the design's decision 3). The tools print the difference; nobody has checked it by hand.
- The flag for open-label extensions (the design's guard against clustering).
- The count deduplicated by lead sponsor. The counts give the number of distinct lead sponsors per class, which is a different number.
- The site.
