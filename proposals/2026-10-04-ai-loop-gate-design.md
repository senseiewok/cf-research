# The AI loop as the gate for review tiers

Status: proposed. Slices 1, 2, 3 and 5 are implemented in the same release train (scripts in the control repo, ledger in this repo); slices 4, 6 and 7 are board tasks. Date: 2026-10-04. Owner repo: hq (the scripts and rules); this note lives here because proposals do.
Provenance: designed by a read-only frontier-model consultant (model override `fable`) from the evidence in `ledger/2026-10-04-first-session.md`, then checked by the controlling agent. Scripts were drafted by the local worker (`qwen3.8:27b-64k`) or written by the controlling agent, and each has an independent test. Parts marked (P) are untested proposals.

## Why

One long working session produced 23 logged errors. Most cost a few steps; one broke a promise to the human (a protected file swept into a commit); the most expensive were claims that went beyond the evidence, and those were caught mainly by a different model reading blind. The loop below is the smallest set of gates that would have caught each class earlier. The evidence is one session by one agent; trust the design after several ledgers, ideally from more than one agent.

## The loop as a gate

One unit of work: scope (tier, exact paths, a done-when command, and for elevated and above a list of factual claims), ground, draft (local worker or orchestrator), deterministic checks, blind model review, orchestrator adjudication (each finding checked against the file), pre-commit gate, human gate, record the acceptance output on the board task.

| Gate row | Routine | Elevated | Full |
| --- | --- | --- | --- |
| Staged set equals the expected set; commit with explicit paths | yes | yes | yes |
| Privacy scan with a canary per pattern | yes | yes | yes |
| Staged scripts and data parse; owning checker passes | yes | yes | yes |
| Counts generated, not typed; rendered-glyph count where a convention names one | when present | yes | yes |
| Local worker review of the diff (findings schema, 2 fast samples) | diff over about 20 lines | yes | yes, never on secret material |
| Claim ledger: type, check, evidence id per claim (P) | no | yes | yes |
| Local worker thinking pass over the claim ledger (P) | no | 1 sample | 1 sample |
| Orchestrator re-reads the complete sentences behind every only, none, every, differs or absent claim | no | yes | yes |
| Blind different-family challenger | no | when public-facing or a factual claim changed (P) | always; unavailable means blocked |
| Human | approves the commit | approves the commit | reads the diff |

`run-gate.ps1 -Tier <tier> -Expected <paths>` runs the deterministic rows and prints the rest as a checklist. It never ticks a model or human row.

Stop rules (unchanged from the council skill): two fast attempts, then one thinking attempt, then the orchestrator does it. Reviews are never retried. One approved cloud call per unit.

## Grounding

| Claim type | The check that settles it |
| --- | --- |
| A quote exists | exact substring in raw text at the cited place; the complete sentence, never a prefix |
| Absence, only, none, every, differs | search the raw text of every unit the claim ranges over, never a parsed structure; the verifier must first find a planted instance |
| Counts in prose | pasted from checker output |
| A convention (the 65 roses) | a script counts rendered instances (`count-rendered.ps1`) |
| Privacy and secrets | `check-staged.ps1`; fails unless each pattern matches its own canary |
| State of git | `git diff --cached --name-only` against the expected set before; `git show --stat` after |
| Anything about models | measured counts with sample sizes; comparative words are dropped without a cited measurement |
| A command did something | its exit code and output |

A verbatim-quote check proves the characters exist. It does not clear an inference: "returns to only the Leeds criteria" is a claim about what the rest of the sentence does not say. A verifier that has never been seen to fire has proved nothing, so each one gets a positive control.

## Diversity and delegation to the local worker

What was measured for the local worker is recorded in the model profile skills (small samples: PowerShell check writing, planted-defect review, an injection probe). Not measured: prose review, inference over registry text, privacy scanning, images. Fast and thinking are one model, so they are not diversity. Gemma 4 31B and Laguna XS 2.1 are installed locally but not admitted as reviewers; they would be benchmark subjects first, and (P) tool-free shadow reviewers after the injection probe, never counted as the challenger.

Goes to the local worker: bounded drafts with a done-when command; findings-schema review of a diff; mapping spans to fixed ids with abstention; a thinking pass asking "which claims does the evidence not support". Does not: deciding done, privacy scans, git state, counts, model comparisons, images, secret material, PHI, deployment specifics, or satisfying the full-tier challenger.

Independence: reviewers see the packet before any findings; the packet carries source, artifact and constraints, never another reviewer's findings or the worker's reasoning; the orchestrator unions and verifies, and never counts votes.

### Delegation results from building the slices (two scripts and two review passes, one session)

| Script | Local attempts | Independent verifier | What the orchestrator then did |
| --- | --- | --- | --- |
| `count-rendered.ps1` | accepted on attempt 2 (fast): 13 of 13 | written first, by the orchestrator | read the code, found an empty-fence edge case, added the case, watched it fail, fixed it (15 of 15) |
| `check-changed.ps1` | 3 attempts (2 fast, 1 thinking): 16 of 17 | positive control: always-pass and always-fail stubs scored 5 of 17 and 1 of 17 | the one failure was my spec asking for `ConvertFrom-Json`, which accepts a trailing comma; switched to `Test-Json`, silenced git stderr, aligned output, added a skip message (17 of 17) |
| ledger critique | one thinking pass | ids checked against the ledger | all four named gaps were real gaps in the first mechanism list |
| diff review of two skill files | two fast samples, same finding in both | each finding checked with a command | the one finding (a measured PowerShell behaviour called wrong) was a false positive in both samples; I re-measured the behaviour (E24) |

An early run wasted three attempts because the orchestrator's code extractor cut the script at its own inner backticks (E20). Two samples agreeing was no evidence on the diff review. The local worker produced working first drafts quickly (6 to 21 seconds per attempt); the defects that mattered were in the spec and the verifier or in edge cases the worker was never asked about. That supports "draft locally, verify independently, read the code", and says nothing about the worker's quality on tasks it was not tried on.

## Improving over time

- **Ledger.** `ledger/YYYY-MM-DD-<topic>.md`, one table per session, written by the controlling agent from its own transcript, read by a human before commit, scanned for private detail. `ledger/tally.py` counts entries by first detector and class.
- **Promotion, one path.** A single observation stays in the ledger. A lesson becomes a guard or a skill line only after a reproduction (a second occurrence, or a fixture that fails without the guard and passes with it), via a board task with `source_evidence` naming the ledger rows, a human setting `ready`, and a `VERSION` bump. The `ConvertFrom-Json` lesson in this session has a fixture, so it qualifies; the console-encoding lesson recurred several times.
- **Expiry.** A guard that has blocked nothing across six months of ledgers is reviewed for removal. A lesson with a `checked` date older than six months is re-verified or marked `[unverified]`.
- **Metrics (P).** Escapes by first detector, rework steps per unit, gate blocks and false blocks. Not worth trusting before several sessions from more than one agent.
- **Handoff.** A new session reads `AGENTS.md`, the `in_progress` board rows, any `Status: proposed` proposals, then the newest ledger.

## Slices

| # | Slice | Status |
| --- | --- | --- |
| 1 | `security-git/scripts/check-staged.ps1` plus `privacy-patterns.txt` (E8, E9, E18) | done; 12-check self-test |
| 2 | `ai-loop-council/scripts/check-changed.ps1` (E1, E2, E11, E21) | done; 17-case test |
| 3 | `ai-loop-council/scripts/count-rendered.ps1` (E16, E17) | done; 15-case test |
| 5 | `ledger/` folder and `tally.py` | done; self-check |
| 4 | claim checker for quote, absence and count claims with a mandatory positive control (E12, E13) | board task |
| 6 | tier table in `AGENTS.md`, slimming `ai-loop-council/SKILL.md` (about 7,300 words), retiring the loop-orchestrator prototype | board task, human gate |
| 7 | (P) admission run for the two installed local models | board task, full tier |

## Risks and unknowns

- One session, one agent, mostly documentation work; code and UI sessions may need different guards.
- The agent writes `-Expected`, so the staged-set check catches drift, not intent; the human still reads the printed list.
- Privacy patterns will flag legitimate text. Count the false blocks or people will bypass the gate.
- The local worker's value on prose is unmeasured. Its elevated thinking pass may cost minutes and find nothing; keep it for two sessions, then read the ledger.
- The blind challenger is the most expensive row and the only one that caught the two largest escapes.

## Decisions for the human

1. Require a blind different-family challenger at the elevated tier for public-facing artifacts and changed factual claims (a change to the tier table in `AGENTS.md`)? Recommended: yes.
2. Keep the ledger in this public repo, sanitized and scanned? Recommended: yes.
3. Delete the loop-orchestrator prototype instead of keeping a stub? Recommended: delete; it never ran, and its config names one model as worker, challenger and council, contradicting the different-family rule.
