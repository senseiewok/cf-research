# Trial atlas: withholding classes that fail the accuracy check

Status: proposed, written 2026-10-11. The policy was chosen by a blind four-family council (`claude-opus-5-5`, `claude-sonnet-5-5`, `fable` and `gemma4:31b-it-q4_K_M`, each answering alone from one packet of counts and design facts). Their answers are model output (tier T0); the decision below is the controlling agent's reading of them, not a vote. Implemented in `tools/trial_atlas/check_atlas.py` (board row T-0130).

## The problem

The gate blocks all publication when any shown class fails its per-class precision or recall test. On a fresh 992-row frozen set (lexicon 0.3.1, the sampled strata never used to tune the rules), 23 of 35 scored classes passed or were too few to judge, and 12 failed. About 28 percent of entries carry at least one failing class and about 25 percent carry only failing classes. Blocking everything wastes 23 sound classes; showing the failing counts publishes numbers the lab's own check says are wrong.

## What the council said

| Question | Opus | Sonnet | Fable | Gemma |
| --- | --- | --- | --- | --- |
| Which option | B (withhold failing classes); D (merge classes) only if declared in writing before scoring and judged on new rows; reject C (badge) | B with D only if pre-declared; reject A and C | B, D as a follow-up under two conditions; reject C | B |
| Show withheld names or counts | names, plus one aggregate of studies only in withheld classes; no per-class counts | same | same, marked approximate | names, no counts |
| The "not shown" share | record and print it next to the counts; a separate cap, set by a person before the next run | compute S2 on the full set and on the shown set | record it and state it on the page; a cap is the maintainer's decision | treat withheld as unclassified |
| Who computes the list | the gate, from the metrics; the page build fails if its list differs | the gate; the page recomputes and fails on a difference | the gate; no hand edits; a person approves by report hash | the gate records a hash |

Risks all four named: the shown page over-represents easy-to-tag classes (patient-reported, microbiology and feasibility outcomes would be missing, and readers may infer trials rarely register them), and a pass is agreement with model labels, approved but not read row by row, with no clinician having read the class definitions. Two warned that a rule round that knows which classes failed must be checked on a fresh sample or it measures tuning.

## Decision

1. **Option B**: withhold, class by class. A shown class that fails precision or recall where estimable, or has no labelled positive, is withheld instead of blocking. A class that is "few checked" stays shown with its note. Every other stop rule is unchanged and still blocks.
2. **Names yes, counts no.** The page lists withheld classes by name and gives one aggregate: how many studies fall only in withheld classes. It does not show a per-class number for a withheld class.
3. **The gate computes the list.** `check_atlas.py --withhold-failing` derives it from the metrics and the thresholds; there is no way to name a class by hand. `--record FILE` writes a counts-only run record with the hashes of the snapshot, lexicon, tags, model tags, key, sheet wording and exclusion list, the thresholds, every class's figures and status, the sorted withheld list and its hash, the not-shown figures and the verdict. `--check-record FILE` verifies a record's internal hashes and that each withheld class really fails the recorded thresholds.
4. **The page build verifies.** The generator reads the withheld list only from a record, calls the check, and fails if what it would show or withhold differs from the record or if the current lexicon hash differs from the record's. (This part lives in the private site repository and is not built here.)
5. **A person approves three things**, by committing the record in a pull request: the thresholds (already set), the reference labels (approved 2026-10-11), and the specific record by its hash. Any change to the list means a new check, a new record and a new approval.
6. **The not-shown share is printed, not capped.** Today 24.9 percent of entries carry only failing classes. A cap is a decision for a person, made before a run and not after seeing the number. The page states the share next to the counts, not in a footnote.
7. **Merging classes (option D) is not done now.** It would need the merge written down before any scoring and judged on rows not used to choose it. This frozen set has been seen.
8. **Page wording for a withheld class** (to be approved by a person): "We counted this class but do not show it: on a check against model-made labels, our tagging fell below the lab's accuracy threshold, so the number is not reliable enough to publish."

## What this does not do

It does not make the shown classes correct. A pass means agreement with a council of models. It does not read any rows, and it does not replace a clinician reading the class definitions. It does not decide the cap on the not-shown share, the page layout, or whether the page is published at all.

## First real run (2026-10-11)

Run on the redrawn census frozen set: lexicon 0.3.2-draft, the consensus-filtered model tags (88 proposed, 70 kept), 938 rows (797 labels reused, 141 new; the sampled strata exclude the 926 development entries, hash checked), the COPD control, 4 planted non-CF ids and 3 canaries. The record is `tools/trial_atlas/run-records/2026-10-11-run-1.json`; `check_atlas.py --check-record` passes 44 of 44 checks on it.

- **Verdict: allowed with 10 withheld classes.** 23 classes shown, 2 few checked (`bone_density`, `not_stated`), 10 withheld.
- **Withheld:** `feasibility_adherence` (precision 0.812), `gi_pancreatic` (0.811), `healthcare_use` (recall 0.691), `liver` (precision 0.818), `lung_function_unspecified` (0.810), `microbiology_sputum` (0.780), `physical_activity` (recall 0.585), `pro_named_other` (recall 0.797, 0.003 below the line), `pro_unnamed` (precision 0.360, recall 0.637), `survival_transplant` (precision 0.800).
- **Not shown:** 388 of 2,052 entries (18.9 percent) carry only withheld classes; 188 studies fall only in withheld classes. No cap is set on this share; that is a person's decision made before a run.
- **S2 (other plus unclassified)** is 5.4 percent against the 15 percent limit.
- **Two things to look at before anyone relies on this list:** `pro_named_other` is withheld by 0.003, and `physical_activity` (12 tagged entries) moved from passing to withheld when the development entries left the sampled strata, so its recall rests on few rows. A different seed could move such borderline classes. The record names the lexicon hash, so any rule change makes it stale.

The reference labels are a council of models' (approved as the reference by the maintainer on 2026-10-11, not read row by row); no clinician has read the class definitions. Nothing is published.
