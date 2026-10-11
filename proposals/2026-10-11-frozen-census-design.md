# Trial atlas: a census-based frozen set for the accuracy gate

Status: proposed, written 2026-10-11. The design was chosen by a blind three-family council (below). The labels it will be scored with are made by models, so every number it produces is tier T2 at most until a person has checked a sample (see "What a person must still do"). Implemented in `tools/trial_atlas/` (board row T-0130).

## The problem

The accuracy gate (`check_atlas.py`, stop rule S1) wants every shown class to have at least one labelled positive in the frozen set, and it measures precision and recall for a class only when 10 or more entries are tagged and 10 or more are positive. The real run has 36 shown classes and 2,052 entries, but the classes are small: 20 of the 36 hold 28 or fewer tagged entries. A random frozen set of 100 to 150 rows cannot give them positives (on 100 rows, 5 classes had none and 1 of 36 classes was estimable), and `make_labelling_sheet.py` caps the set at 150 rows. The result is a gate that cannot pass, whatever the tools' real accuracy.

## What the council was asked and said

Question packet: counts only (no registry text). Members: `claude-opus-5-5`, `claude-sonnet-5-5` and `gemma4:31b-it-q4_K_M` (three families; each answered alone). Their answers are model output, kept as leads; the decision below is the controlling agent's reading of them, not a vote.

| Point | Opus | Sonnet | Gemma | Decision |
| --- | --- | --- | --- | --- |
| Design | (b) census of small classes, plus a sample of large ones; reject (c) | (b) plus (a)'s targeted extension; reject (c) | (b) | (b), and keep S1's "no labelled positive" block |
| Where the census stops (N) | 30 | 30 | 50 | **30**: in the real counts the class sizes run 28, then 32 (`not_stated`, which has its own census stratum), 40, 47, 54, so 30 takes exactly the 20 classes at 28 or fewer. Gemma's 50 would add classes of 40 and 47 entries, whose censuses cost more and add little over a 40-row sample. Two of three members said 30 |
| Large classes | at least 30 tagged rows each | 40 to 60 | a random sample | **40** tagged rows per large class |
| Recall | draw from untagged entries and census the unclassified pool; weight by inverse probability; the tagged rows alone give precision only | same, report the interval and the count behind it | sample entries not tagged with the class | a census of the unclassified pool plus a weighted random sample of the entries tagged only with large classes |
| Cap of 150 | remove, the census sets the size | remove | not stated | removed for the census design; the old design keeps it |
| `other` | not scored (the tools never predict it) | residual, not scored | residual | not scored |
| Consensus rule | acceptable only if frozen now and scored on rows labelled after the choice; report the plain route beside it | post-hoc; report plain and consensus side by side, the date of the choice, and the number dropped | report the trade-off | frozen from 2026-10-11 as a rule chosen after the results, with that origin stated (below); the approvals council (`2026-10-11-council-approvals.md`) asked that it not be called pre-declared |

Risks the council named, all three of them in some form: the labels are model labels, so a class's "accuracy" is agreement among models until a person checks a sample; correlated errors of the council and the tools inflate both precision and recall; the frame is defined by the tools' tags, so a change to the rules or the model route makes a census stale and it has to be redone; rare positives in the large untagged pool may never appear in an affordable sample, so a small class's recall stays loosely bounded.

## The design

**Frame.** All in-scope entries (an entry with no primary-outcome wording is never drawn). The classes an entry carries are its final classes: its rule tags plus its accepted model tag, if any. The frame is bound to the snapshot, the lexicon, the rule tags and the model tags by their hashes.

**Strata.** Each entry is in exactly one stratum, assigned in this order:

| Stratum | Which entries | How many are drawn | Weight |
| --- | --- | --- | --- |
| `census:small` | carries any class that has `small_max` (default 30) or fewer entries | all | 1 |
| `census:unclassified` | no accepted class (status `unclassified`) | all | 1 |
| `census:not_stated` | status `not_stated` | all | 1 |
| `sample:<class>` | carries only classes with more than `small_max` entries; the stratum is the entry's rarest such class (ties by the lexicon's class order) | `large_per_class` (default 40) per stratum, or all if fewer; a seeded random draw | population / drawn |

The sheet is blind exactly as before: wording only, no tag, no stratum, rows shuffled. The key is sealed.

**Key** (`design: "census-v1"`; the old key has no `design` field and keeps its meaning): the fields of the old key, plus `design`, `small_max`, `large_per_class`, `model_tags_sha256` (or null), `consensus` (null, or the rule's record, below) and `strata`, an object keyed by the stratum names above, each `{"kind": "census" | "sample", "population": n, "drawn": n, "weight": population/drawn}`. Each row's `stratum` is one of those keys. The sum of the strata's populations is the frame's size, and the gate checks it.

**Estimates** (in `check_atlas.py`, when the key's design is `census-v1`). For a class c, each labelled row counts with its stratum's weight: true positive when c is in both the entry's final classes and the council's labels, false positive when only the final classes name it, false negative when only the labels do. Precision is weighted TP over weighted TP plus FP; recall is weighted TP over weighted TP plus FN. Weights are 1 in the census strata, so a class whose entries are all in census strata has exact precision. Intervals come from a seeded stratified bootstrap (2,000 resamples, seed recorded in the output); the Wilson interval stays for the old design. A class is estimable for precision when 10 or more entries in the frame carry it, and for recall when its weighted positives are 10 or more. The thresholds do not change (precision 0.85, recall 0.80, 10 or more). A class that is not estimable is noted "few checked", as now. `other` is not scored. S1's other rules stay: a minimum of 50 labelled rows, no unlabelled row, no shown class with no labelled positive.

**What the recall number means.** Recall for a small class is exact over the census strata and a weighted estimate over the sampled ones. A positive the tools missed can sit in a sampled stratum without being drawn; the interval says so, and the page's "few checked" note stays where the interval is wide. A recall measured only on rows drawn because a tool tagged them would read near 1 and is never reported as recall.

## The consensus rule, frozen now

A model proposal for an unclassified entry is kept only when a blind model of a different family, labelling the entry without seeing the proposal, also gave that class. It was chosen after the plain route's precision of 0.759 was seen on the development census (0.885 under the rule, which dropped 21 of 108 proposals, 16 of them wrong and 5 right), so that number is a development figure. The rule is frozen with this proposal; the gate prints the plain route and the rule side by side with the count dropped; the page says in one sentence that the rule was chosen after seeing the plain route. Entries labelled after this decision and scored under it are the only confirmation.

## Frozen after labelling

The snapshot, the lexicon version, the model profiles and prompts, the consensus rule, the frame and strata, `small_max`, `large_per_class`, the seed, the thresholds and the class list. A change to any of them is a new frozen set, not an edit. A label changed after the council has answered is a logged decision on a disagreement, never a silent edit.

## What a person must still do

A model council cannot make a number T3. Before anything is published a person approves: the class definitions and labelling guide the council used; a spot-check of the labels (every census row on which the council members disagreed, and a random sample of the rest); the thresholds and their "not a standard" wording; the consensus rule's post-hoc status; the S4 control results; and the page wording, including the statement that the labels come from a model council.

## Not covered here

The non-CF control snapshot (stop rule S4) and the canaries are a separate step that needs a second fetch. The council's `other` label had no tool counterpart and stays unscored. Whether the council's labelling guide matches what a clinician reads in the classes has not been checked by a clinician.
