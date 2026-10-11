---
name: trial-atlas-pipeline
description: The rules for the trial endpoint atlas, the count of what interventional cystic fibrosis trials register as their primary outcome in one dated ClinicalTrials.gov snapshot. Use before any work that touches the atlas data, the snapshot or its fetch, scope, the rule lexicon or its tags, the model-assisted tagging route or its challenger, the frozen labelling set, counts.json, the gate and its stop rules, or anything that publishes or describes the counts. Covers the order of work, the approval gates, who may do each step, where every output lives, the stop rules by name, what a verified quote proves, and how lessons from runs become rules.
license: CC0-1.0
compatibility: Works with any Agent-Skills-spec-compatible tool that reads .claude/skills (Claude Code, GitHub Copilot; see AGENTS.md in cf-lab). The tools it points to need Python 3.10+; the model-assisted route also needs PowerShell 7 and a local Ollama model. A local worker model does not load skills; it gets pasted packet text (see source-grounding/references/local-worker-packets.md).
---

# Trial endpoint atlas pipeline

The atlas counts what interventional CF trials registered as their primary outcome, from one dated, hashed snapshot of ClinicalTrials.gov. The design is `proposals/2026-10-10-trial-endpoint-atlas.md` (board row T-0130). Every command is in `tools/trial_atlas/README.md`; this skill does not copy them. It holds the general rules. Specifics seen on real runs are in `references/lessons.md` until they are promoted. For any sentence about the atlas, load `source-grounding` too.

## How each agent gets this skill

| Agent | How |
| --- | --- |
| Claude Code | Start in `cf-lab` and add this repository with `--add-dir ../cf-research`, or `/add-dir` during a session; AGENTS.md in cf-lab says skills from an added folder reload live. `permissions.additionalDirectories` alone gives files, not skills |
| GitHub Copilot | AGENTS.md in cf-lab says Copilot reads `.claude/skills/`. Whether it reads that folder in every folder of the multi-root workspace was not tested here: see AGENTS.md |
| A local worker model | Loads no skill. The controlling agent pastes a packet from `../source-grounding/references/local-worker-packets.md` and the tools' own templates |

## The order of work

How to run each step: the README section named in the last column.

| Step | Needs | Produces, and where |
| --- | --- | --- |
| 1. Fetch the CF snapshot | A person's approval for the network call; `EVIDENCE_CONTACT` from the environment only, never on a command line or in a file | Pages and a hashed manifest, in git-ignored `sources/downloads/` or the private files folder. A stopped fetch is never retried in a loop. README: "The order of work", "When a fetch stops" |
| 2. Fetch the non-CF control; confirm both | The same approval; a person compares the pages with the registry's documentation before confirming | `config_verified` in each manifest. README: "To confirm on the first real run" |
| 3. Scope and rule tags | An agent or a second person, never the labeller | A private tags file. README: "The order of work" |
| 4. Draw the frozen set | The same runner; a seed chosen and written down at the draw | A blind sheet, its terms file and a sealed key kept apart. README: "The order of work" |
| 5. Label | A person other than the runner of steps 3 and 4, who never sees their output or the key | The labelled sheet |
| 6. Model tags for the unclassified rest | A local proposer and a challenger from a different family, each named explicitly | Private outputs outside every git working tree. README: "The model-assisted route" (arrives with the model-route change of T-0130) |
| 7. Counts | All of the above | A private counts file carrying the registry's terms block. README: "What counts.json holds" |
| 8. The gate | Every gate input, real controls included | Exit 0 only when integrity and every stop rule pass. README: "Exit codes and stop rules" |
| 9. Pages | A separate private step, outside this repository | It builds pages from a counts file and refuses synthetic data |

Class names, glosses and model "family" leads wait for a clinician reader and a council before any class or rule changes. Glosses are aids for the labeller, not public text.

## Privacy and the registry's terms

- Nothing registry-derived is committed: no pages, tags, sheets, keys, model replies or real counts. The catalog entry `clinicaltrials-gov` keeps `redistribute` false until the maintainer decides.
- The registry's terms items travel with every count (the terms block in the counts file, the gate's printout, the sheet's terms file). Never strip or reword them. README: "Scope and the registry's terms".
- No contact address, machine detail or private path in any committed file.

## Stop rules

Named as the README and `check_atlas.py` name them; exit 3 means one tripped. Thresholds can change before the frozen set is labelled, never after.

| Rule | Guards against |
| --- | --- |
| S0 | Synthetic data reaching publication |
| S1 | A missing, unbound (key hashes), incomplete or too-small frozen set, and classes below the precision and recall floors |
| S2 | Too many entries left "other" or unclassified (above 15 percent) |
| S3 | A missing recall route, or an unexplained difference between the two retrieval routes |
| S4 | Controls not run, run only on synthetic fixtures, or too small |

## What a verified quote proves

Only that the words exist in the entry. It does not prove the class is right; a claims-checker pass does not prove a quote supports its sentence; agreement between two models is not correctness. The person-labelled frozen set and the blind challenger measure correctness.

## Rules for agents

- Write the test first and show it fail.
- Read the files a subagent or model says it changed before accepting its claim. No completion claim without the command output that proves it.
- Never put a registry quote in a public file beyond the lab's quotation limit (8 words; `tools/sources/check_source_overlap.py --quotes-only`).
- Never invent a contact, a label, or a seed after the fact. Never run steps 3 and 4 for someone who will label.
- Everything public passes the private page step's prose ledger and registry-claims check, and `check_registry_claims.py` here (README: "Grounding the registry statements").
- Observed fact, computed value and inference stay labelled apart. An inference is never page text.

## Lessons: how specifics become rules

This file is the catch-all. Specifics live in `references/lessons.md` as dated observations (date, what was seen, tier, a "reproduce with:" line, status: observed, proposed-rule or promoted). An observation becomes a rule only through a `proposed` board task, a person setting it `ready`, and a narrow edit to the owning skill with a check (AGENTS.md, "Local files and memory" and "AI loop and delegation"). Delete the observation once promoted. Review or remove an unpromoted observation older than 30 days. Never promote a model's opinion into policy.

## Verification

- The tests: README "Run the tests" (synthetic fixtures, no network).
- The route's planted controls: the self-check in "The model-assisted route".
- A counts file: the self-consistency check in "What counts.json holds".
