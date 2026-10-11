# Trial atlas: council decisions on six approvals

Status: proposed, written 2026-10-11. The maintainer asked for these approvals to go through a council that includes Fable. Every answer below is model output (tier T0 to T2). No model answer is a person's check (T3), and where this note says "approved", it means the council approved and the maintainer's merge of the pull request that carries it is the final step.

Members: `claude-opus-5-5`, `claude-sonnet-5-5`, `fable` and `gemma4:31b-it-q4_K_M` (four families, each answering alone, from one packet of counts and design facts; no registry text). The decisions are the controlling agent's reading of the four answers, not a vote.

## Decisions

| # | Question | Council | Decision | Applied now? |
| --- | --- | --- | --- | --- |
| 1 | The class definitions and rules the labellers used | Gemma: approve. Opus, Sonnet, Fable: approve with changes | Approved as the guide the existing labels were made with. The changes below go into guide and lexicon 0.3.1 for the next draw. They are not applied to the guide now, because an edit after labelling changes what the labels mean. The public page still withholds plain-language class notes until a clinician reads them | Recorded here |
| 2 | The accuracy thresholds (precision 0.85, recall 0.80, 10 or more, at least 50 rows) | All four: approve | Approved unchanged. Wording changes: say the thresholds were chosen by the lab before labelling, are not a standard, and are measured against labels from a council of models, not clinicians. Show each class's labelled count beside its verdict (at 10 entries one error decides the line) | Wording queued for the site generator |
| 3 | The consensus rule and its page sentence | Opus, Sonnet, Fable: not "pre-declared", approve as "chosen after the results, frozen from 2026-10-11". Gemma: approve | Approved as a post-hoc rule, frozen from 2026-10-11, to be confirmed on entries labelled after the choice. Final sentence below | Proposal wording fixed; page sentence queued |
| 4 | A second ClinicalTrials.gov fetch for the non-CF control | All four: approve with changes; all four chose COPD | Approved, bounded as below | Not run yet |
| 5 | The lab's name on the atlas pages | All four: "Sensei Ewok's CF Lab" | Recommended, not applied: it changes the site-wide header brand, and three of the four also advise replacing the "Cystic Fibrosis AI Lab" subtitle that the maintainer chose. Both are the maintainer's calls, so they go in a separate site change | No |
| 6 | Show percentages on the class pages | All four: show, with limits. Opus: drop the percentage card, keep the column | Approved with the limits below | Queued for the site generator |

## Item 1: changes queued for guide and lexicon 0.3.1

Raised by at least two of the three reviewers who asked for changes, unless marked.

- `inflammation_markers`: remove "exhaled nitric oxide is mainly a type-2 marker and usually low in CF". It is a clinical claim, not a labelling rule, and has no source here (all three).
- `liver` against `safety_tolerability`: "liver blood tests are safety" mislabels a trial whose endpoint is a liver enzyme in liver disease. Say a liver blood test named as the liver-disease endpoint is `liver`, and a routine safety lab is `safety_tolerability` (Opus, Fable; Sonnet asks for the same pointer for glucose and bone).
- `exacerbations` against `healthcare_use`: an entry that calls it an exacerbation is `exacerbations`; one that counts courses or admissions without that word is `healthcare_use`. Drop "in CF" from the gloss so it holds on a non-CF snapshot (Opus, Fable).
- `oxygenation`: end-tidal carbon dioxide is ventilation, not oxygenation; rename the class or drop it from the gloss (Opus, Fable). Say exhaled gases for oxygen or carbon dioxide are `oxygenation` and exhaled nitric oxide is `inflammation_markers` (Sonnet).
- `spirometry_other` against `lung_function_unspecified`: the gloss says "or spirometry with no FEV1 named" while the rules say a vague "lung function" is unspecified; say the FEV1/FVC ratio is `spirometry_other`, and that "ppFEV1" and "percent predicted FEV1" count as FEV1 named (Sonnet).
- A tie rule for a measure that fits a specific and a general class: give the specific class only. Questionnaire-based measures against `pro_named_other`: label both only when the entry names both (Opus, Sonnet).
- Say several independent measures in one entry are not a composite (Sonnet).

## Item 3: the consensus rule

It was chosen after the plain route's precision on the development census was seen (0.759; 0.885 under the rule), so it is not a rule declared in advance. It is frozen from 2026-10-11. Page sentence, edited from the four answers:

"A model's suggestion for an entry the rules could not classify is counted only when a second model from a different family, labelling the entry without seeing the suggestion, gave the same class. We chose this rule after seeing the first model's results on the development entries, so its accuracy on new entries has not yet been confirmed. It also discards some correct suggestions (5 of 108 in development), and two models can share the same mistake."

## Item 4: the control fetch

- One non-CF condition: COPD. Bronchiectasis can include or mirror CF trials; asthma skews to children and exhaled nitric oxide.
- Fetch only what the gate needs. It needs 50 outcome entries, so at most 5 requests at 100 studies per page, not 20.
- Leave out any study that also lists cystic fibrosis as a condition.
- Same client, rules and terms re-read as the first fetch; public API v2, honest user agent, a pause between requests, no retry past a refusal, robots.txt checked, HTTPS and same-host redirects only; stored in a private folder outside every repository, never redistributed.
- Record the snapshot date and the query.
- Write down before the fetch which classes the control watches. They are the gate's own control set, the CF-specific classes the lexicon marks (sweat chloride, nasal potential difference, intestinal current, other direct CFTR measures), plus any class the lexicon marks as a negative control. Classes common in COPD (FEV1, exacerbations, healthcare use, patient-reported outcomes) do not count against the gate.

## Item 6: the percentage column

Show the column; drop the percentage card (a headline percentage invites reading it as a success rate). Counts first. "k of n" under 10 studies, as now. For a class that is "few checked", counts only. Label the column "share of that year's studies registering it" and say what the denominator is. Flag recent start years as incomplete in this snapshot. Keep the note that a larger share does not make a measure better.

## A fourth-family audit of the existing labels

Fable had labelled nothing before. It labelled 153 of the 427 existing council-labelled entries blind: all 73 on which the council's models did not agree, and 80 drawn at random (seed 20261011) from the 354 on which they did. Same instructions and class list as the council; it saw no council label. Counts only:

| Audited rows | n | Fable's labels equal the council's set | Share at least one class |
| --- | --- | --- | --- |
| Rows the council agreed on (random draw) | 80 | 77 (0.96; 95% CI 0.895 to 0.987) | 78 |
| Rows the council's models split on | 73 | 55 (0.75; CI 0.644 to 0.838) | 61 |
| Both | 153 | 132 (0.86) | 139 |

The disagreements are mostly one pattern: the council's `other` against a specific class from Fable (`nutrition_growth` 5 times, `inflammation_markers` 2). That fits the council's `other` label being the least stable one, which is also why `other` is not scored. This is one model's agreement with the council on a sample, so it is tier T2 for the sample and says nothing about the 274 rows not audited; it is not a person's check.

## What a person must still do

**Labels: approved by the maintainer on 2026-10-11, as the reference for the gate.** The approval relies on the council runs (Opus, Sonnet and Gemma labelling, Fable and a model adjudicator on the splits) and on the fourth-family audit above. It is an approval of the method and its record, not a row-by-row check: no person has read the individual rows, so the numbers built on these labels are still labelled as measured against model-council labels. A spot-check from the private review pack stays open for anyone who wants one. Have a clinician read the class definitions before any plain-language note is published. Decide the lab's name and the header subtitle. Merge.
