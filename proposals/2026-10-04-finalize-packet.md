# Task packet: finalise the 2026-10-04 changes and prepare PRs

Status: applied on 2026-10-04 in local commits (hq `f1fe8ff`, research `d0cde10`), not pushed or merged. Board tasks T-0003, T-0004 and T-0022 stay `in_progress` until a human sets `done`. Original text follows unchanged.

Paste everything below the line into the orchestrating agent in `hq`. It follows the `ai-loop-council` packet contract and supersedes `2026-10-04-implementation-packet.md`, which it includes. Review tier is **Elevated** throughout (agent instructions are being edited) and **Full** for the license step.

---

## Goal

Four repos have uncommitted work from an external research session dated 2026-10-04. Bring each to a reviewable pull request: apply the skill-update proposal to `hq/.claude/skills/`, adopt the task board, integrate the new research documents, and leave nothing in a half-applied state. Do not extend scope beyond what the proposal files and board tasks name.

## Grounding: read in this order, do not skip

1. `hq/AGENTS.md`, including the new "Task board" section.
2. `hq/tasks/README.md`, then `hq/tasks/BOARD.md`. The tasks you are executing are **T-0003, T-0004, T-0022**, plus the PR preparation below. Everything else on the board stays `proposed`; you do not start it.
3. `../research/proposals/2026-10-04-hq-skill-updates.md`: the exact skill text to apply, item by item, with evidence.
4. `../research/proposals/2026-10-04-implementation-packet.md`: the step list and acceptance checks for the skill edits. Apply it as written.
5. `../research/sources/README.md` and `../research/sources/catalog.yaml`: what the edited skills will point at.
6. `../research/proposals/2026-10-04-website-review.md` and `2026-10-04-youtube-about.md`: context only; their tasks are not in this packet.

## Constraints

- **Apply, do not author.** Skill text comes verbatim from the proposal. If a passage conflicts with the current skill because the skill changed, stop and report; do not reconcile it yourself.
- **No new facts, no fetching.** Every inserted claim already carries a source and a check date. The only network call permitted is `fetch_sources.py` in the smoke test, which reaches `pr.ecfs.eu` only.
- **Never touch `www.cff.org`, `databases.lovd.nl`, or `cftr2.org` programmatically.** Not with any client. Reasons are in `research/landscape/registry-data-access.md` and are not re-evaluated here.
- **A human approves every diff under `hq/.claude/`.** Stage it, show it, wait.
- **Frontmatter stays valid**; `name` matches folder; `license` and `compatibility` unchanged.
- **Do not push. Do not merge.** Branches and PR descriptions only; the human opens and merges.
- **Do not edit `board.json` statuses to `done` yourself.** Set `in_progress` when you start a task; a human sets `done` after seeing the acceptance output.
- **No PHI, secrets, machine details, or local paths** in any file or PR body.
- **The board is the only task list.** If you find tasks you think are missing, append them to `board.json` as `proposed` with `created_by.kind: model` and your `model_id`, and say so in the report. Do not create a second list anywhere.

## Steps

### A. `hq` repo

1. Set T-0003, T-0004, T-0022 to `in_progress` in `tasks/board.json`; re-render.
2. Create branch `chore/2026-10-04-registry-access-skills`.
3. Apply the four skill edits exactly per the implementation packet steps 1-4 (`security-browsing`, `cf-research-context`, `ai-loop-council`, `model-onboarding`).
4. Gate `.github/loop-orchestrator/loop.py` with a hard failure per implementation packet step 6; update its mention in `AGENTS.md`.
5. Add `tasks/render-board.ps1 -Check` to whatever pre-commit or CI mechanism the repo already uses. If none exists, add a one-line note under "Validation expectations" in `.github/copilot-instructions.md` and do not invent a CI system.
6. Add `cf-projects/CF-Project-Ideas.md` scope-triage table review to the PR description as a **human review item** (T-0020); do not edit the table.
7. Stage everything. Produce `git diff --stat`. Run acceptance checks (below). Write the PR description to `hq/.github/PR-2026-10-04.md` using the template at the end of this packet.

### B. `research` repo

1. Create branch `feat/source-catalog-and-proposals`.
2. Run the smoke test from the implementation packet step 5. Confirm `git status --porcelain` shows nothing under `sources/downloads/`.
3. Confirm `README.md` layout table lists `sources/`, `tools/sources/`, `proposals/`, and `landscape/registry-data-access.md`.
4. Stage everything. Write the PR description to `research/PR-2026-10-04.md`.

### C. Website repo

No code changes in this packet. Create `PR-NOTES-2026-10-04.md` at the repo root listing board tasks T-0014 to T-0018 and T-0023 with one line each, so the website owner sees them in the repo and not only on the board. Nothing else.

### D. `skills` repo

No changes. Record in the report that it was checked and untouched.

## Acceptance checks (verifier decides, not a model)

Paste real command output for each. A sentence saying it passed is a failure of this packet.

- [ ] `pwsh -File hq/tasks/render-board.ps1 -Check` exits 0 and prints `BOARD.md is current`.
- [ ] `python -c "import json; json.load(open('hq/tasks/board.json'))"` exits 0.
- [ ] `python -c "import yaml; yaml.safe_load(open('research/sources/catalog.yaml'))"` exits 0.
- [ ] A script asserts, for each of the four edited `SKILL.md` files: file starts with `---`, frontmatter parses, `name` equals folder name. Output shows four passes.
- [ ] `python research/tools/sources/fetch_sources.py --manual` lists only `access: manual` catalog entries (the Cystic Fibrosis Foundation reports and CFTR2) and no `pr.ecfs.eu` URL.
- [ ] `python research/tools/sources/fetch_sources.py --only ecfspr-adr-2022` then `--verify` prints `1 match, 0 changed, 0 missing`, exit 0.
- [ ] `git -C research status --porcelain` contains no `sources/downloads/` and no `*.part`.
- [ ] `python hq/.github/loop-orchestrator/loop.py` exits non-zero.
- [ ] `git -C hq diff --stat main` touches only: four `SKILL.md` files, `AGENTS.md`, `.github/loop-orchestrator/loop.py`, `.github/copilot-instructions.md` (if step A5 used it), `tasks/*`, `cf-projects/CF-Project-Ideas.md` (already modified; carry it), `.github/PR-2026-10-04.md`. Anything else is a scope violation; unstage and report it.
- [ ] A grep for the old drive-letter machine path across `hq`, `research` and the website repo returns nothing in the files this packet touched.
- [ ] No file in any staged diff contains an email address, an API key pattern, or an absolute `C:\Users\` path.

## Escalate to the human, and stop, if

- A proposal passage conflicts with current skill text (skill changed since 2026-10-04).
- Any acceptance check fails twice after one bounded repair.
- Any step would need network access beyond `pr.ecfs.eu` or a package outside `research/tools/sources/requirements.txt`.
- You are tempted to improve wording in a skill, a proposal, or the website. That is a separate task; append it to the board instead.

## Report back with

1. The three `git diff --stat` outputs (`hq`, `research`, and the website repo).
2. Command output for every acceptance check, in order.
3. The paths of the three PR description files.
4. One line per board task touched: id, new status, and the acceptance output that justifies it.
5. One line per anything skipped or appended to the board, with the reason.

No narrative summary of what the skills now say. The human reads the diff.

---

## PR description template

```
## What

<one paragraph: which board tasks (ids) this PR completes and what it changes>

## Why

<two or three sentences; link the proposal file(s) under research/proposals/ this applies>

## Evidence

<bulleted list: each acceptance check and its real output, abbreviated to the deciding line>

## Human review items

<bulleted list of things a model should not decide: e.g. T-0020 scope triage, any skill wording the reviewer should read in full>

## Provenance

Proposed by: <created_by from the board tasks, model_id or human name>
Applied by: <this agent's model_id>, <date>
Review tier: <routine | elevated | full>

## Not in this PR

<board task ids deliberately left proposed, one line each>
```
