# Document organization and memory design for the lab repos

Status: proposed. Tier 0 and Tier 1 below are done (2026-10-04). Tiers 2 to 5 wait for human decisions and are board tasks T-0029 to T-0032, all `proposed`.
Date: 2026-10-04. Owner repo: hq (the rules); this note lives here because proposals do.
Provenance: drafted by a read-only frontier-model consultant (model override `fable`), then spot-checked by `claude-sonnet-5-5`. Findings marked verified were re-checked against the files; findings marked unverified were not.

## Why

There are 109 markdown files across the lab repos. Most are not prose to maintain: 53 are frozen test fixtures under the council skill. The real problems are small and specific: missing or stale status lines, a rule document that contradicts itself, two always-loaded instruction files that can drift, and committed PR descriptions that go stale at merge. This design fixes those with the fewest new mechanisms.

## Part 1: organization rules

| Type | Lives | Header it needs | Lifecycle |
| --- | --- | --- | --- |
| Always-loaded instruction | `hq/AGENTS.md`; `hq/.github/copilot-instructions.md` | version in `VERSION` | elevated review; each change bumps `VERSION` |
| Skill | `<repo>/.claude/skills/<name>/SKILL.md` | `name`, `description` (existing rule) | draft, reviewed, superseded |
| Test fixture | `ai-loop-council/cases/**` | frozen once its reference check passes | leave as is |
| Task | `hq/tasks/board.json` (`BOARD.md` is generated) | full schema, enforced by the renderer | proposed, ready, in_progress, done, blocked, dropped |
| Proposal, evidence, provenance, idea note | `research/proposals/`, `research/landscape/`, research root | **first line `Status:`** with one of `proposed`, `applied (repo, commit or task id)`, `superseded by <file>`, `review note`, `provenance note`, `archived`; plus a date | proposals are kept after they are applied; git history is the record |
| PR description | nowhere durable | none | delete at merge; keep provenance in board task `notes` and git history |
| Generated file | `tasks/BOARD.md`, the workspace file | "do not edit" header or git-ignored | regenerate |

No new tooling is proposed. The one enforceable check is `grep -L '^Status:'` over the proposal folder.

### Problems found (verified unless marked)

| Path | Problem |
| --- | --- |
| `hq/.github/PR-2026-10-04.md`, `research/PR-2026-10-04.md` | PR descriptions committed in repos; stale at merge |
| `research/proposals/2026-10-04-hq-skill-updates.md` | said "not applied" after it was applied. **Fixed in Tier 1** |
| `research/proposals/2026-10-04-implementation-packet.md`, `finalize-packet.md` | no `Status:` line. **Fixed in Tier 1** |
| the provenance note in `research/proposals/` | no H1. **Fixed in Tier 1** |
| `hq/README.md` structure tree | omits `tasks/` and `Modelfile` |
| `hq/.github/loop-orchestrator/README.md` | advertises "consensus voting", contradicting `AGENTS.md` ("never decide a finding by majority vote"); the prototype is gated but its README still reads as a working feature |
| `hq/.github/copilot-instructions.md` | overlaps `AGENTS.md` (sibling table, skills rule, routing); two always-loaded sources can drift |
| `hq/.claude/skills/ai-provider-compatible-skills/SKILL.md` | duplicate section numbers (7 and 8 twice), a `.skills/` directory tree that contradicts its own location rule, and a generic "AI loops" section |
| `hq/cf-projects/CF-Project-Ideas.md` and `research/cf-ai-ideas.md` | two idea inventories in two repos. A cross-repo move needs a human decision |
| `hq/tasks/` tasks that describe the website | the website repo's board key is kept; board text carries no deploy details or paths (decided 2026-10-04) |
| `skills` repo | no commits, two untracked files; whether that is intentional is unverified |

Leave as is: the skill layout and fixtures, `board.json` and `BOARD.md`, the source catalog and its README, dated landscape notes, `SECURITY.md`, `LICENSE.md`.

## Part 2: memory design

| Layer | Holds | Budget | Must not hold |
| --- | --- | --- | --- |
| L0 `AGENTS.md` | rules, catalog, review tiers, pointers | cap about 2,500 words (now about 2,000) | procedures, evidence, model scores |
| L0b `copilot-instructions.md` | Copilot-only items and "read AGENTS.md" | cap 40 lines | anything also in `AGENTS.md` |
| L1 skills | procedures, profiles, lessons with a verification step | 500 lines each | private paths, transcripts, session ids |
| L2 board | durable intent, handoff state (`in_progress` plus `notes`), measured outcomes | one task per intent | packets, drafts |
| L3 evidence and decisions | dated, labelled facts; applied and superseded proposals | `Status:` line mandatory | undated claims |
| L4 per-user agent memory (outside the repos) | machine-specific operational facts only: where the repos are, which local profile is in use, processes not to stop | a few short files | lab rules or lessons, task state, secrets, PHI, claims about which model is running |
| L5 local model | nothing; it gets a fresh packet built from L0 to L3 | one packet | earlier packets |
| Ephemeral | PR descriptions, `attempts.jsonl`, logs, transcripts | git-ignored or deleted | n/a |

**Capture and promotion, one path.** Observation, reproduce it with a cheap check, board task `proposed` with `source_evidence`, human sets `ready`, fix, the lesson lands narrowly in the owning skill with a fixture or check where the mistake recurs, then a `VERSION` bump. An observation that cannot be reproduced stays in the task's `notes` and is never promoted. Models propose; they never write L0 or L1 directly.

**Expiry.** A `checked <date>` older than six months is re-verified or relabelled `[unverified]`. Applied and superseded proposals stay. PR descriptions are deleted at merge. L4 entries are deleted when the path or process changes.

**Handoff.** A new session reads `AGENTS.md`, then the `in_progress` rows of `BOARD.md` with their `notes`, then any `Status: proposed` file in `research/proposals/`.

**Never stored in shared memory:** secrets, PHI, drive-letter paths, email addresses, session ids, raw model output, third-party instructions, private site domain or deploy details, a model's statement about its own identity.

**Human audit** (quarterly, or at each `VERSION` bump): `git log -p AGENTS.md .claude/`; `render-board.ps1 -Check`; `grep -L '^Status:' research/proposals/*.md` must be empty; grep the public repos for the site domain, drive-letter paths and `@`; read the L4 folder end to end.

## Part 3: tiers

| Tier | Actions | Repo | Executor | Review | Acceptance check | Not done |
| --- | --- | --- | --- | --- | --- | --- |
| 0 measure | inventory, link and status greps, privacy grep | all | script, frontier consult | n/a | exit codes recorded | any edit |
| 1 status lines | add or fix `Status:` on all six proposals; one H1 each | research | deterministic edit; no model needed for six headers | routine | `grep -L '^Status:'` empty; exactly one H1 per file | moving or renaming files |
| 2 hq housekeeping | README tree adds `tasks/` and `Modelfile`; first line of the loop README says "not functional, see AGENTS.md"; after merge delete the PR description and copy its provenance to task `notes`. **Human gate: `.github/` is protected** | hq | local worker drafts; human commits | routine | `grep -c 'tasks/' README.md` is 1 or more; the first five lines of the loop README contain "not functional"; `render-board.ps1 -Check` exit 0; no `PR-*.md` tracked | `loop.py`, settings, skills |
| 3 instruction dedup | shrink `copilot-instructions.md` to a pointer plus Copilot-only items; fix the authoring skill (unique numbering, drop the `.skills/` tree and the AI-loops section). **Human gate: agent instructions** | hq | cloud or frontier drafts; human approves the diff | elevated | file at most 40 lines; unique numbering; `check-skill-frontmatter.py --against main` passes; the setup self-test passes | changing what any rule means |
| 4 cross-repo | move the idea inventory to `research`, one commit per repo; first commit in `skills`; decide the private-repo naming. **Human only** | hq, research, skills | human | elevated | file exists in exactly one repo; links resolve | publishing any skill |
| 5 memory hygiene | tidy per-user memory (outside the repos) | none | human | n/a | folder listing only | writing lessons into L4 |

Tiers 1 and 2 are independent. Tier 3 waits for the current hq branch to merge, so `--against main` is meaningful.

## Decisions for the human

Answered 2026-10-04: all three "yes" as recommended (decision 3 after the T-0020 review).

1. Make `AGENTS.md` the only rule source and reduce `copilot-instructions.md` to a pointer? Recommended: yes.
2. Stop committing PR descriptions and keep provenance in board fields and `Status:` lines? Recommended: yes, delete the files at merge.
3. Move the idea inventory to `research` after the human review of its triage table (T-0020), and keep the website repo's board key while writing "no site domain, paths or deploy details in board text" into `tasks/README.md`? Recommended: yes.

## Unverified

Whether the empty history in `skills` is intentional; the contents of per-user memory (not read by the author of this note).

## Update 2026-10-04: where layer L4 lives (applied in hq)

A read-only frontier design review (no files written, no network) fed this update. What was applied, in hq `AGENTS.md` ("Local files and memory") and `setup.cmd` / `setup.sh`:

- **L4 is now `cf-lab-files/memory/`**, a folder beside `hq` (now `cf-lab`), not a git repo, never published, created by the setup scripts with a README that lists what must never go in it. It is the only lab-level untracked memory, so Copilot and a Qwen worker can read it as well as Claude Code. Each tool's own per-user memory keeps one pointer line and tool-specific preferences; it is not mirrored, because two copies drift.
- **Contents:** `handoff.md` (current session notes, overwritten), `machine.md` (paths, local model, processes not to stop), `observations.md` (dated candidates, each with a "reproduce with:" line). Files are created on first use; nothing else is created, and `downloads/` was dropped (the fetcher keeps `research/sources/downloads`).
- **Reading order for a new session:** `AGENTS.md`, the `in_progress` rows of `BOARD.md` with their notes, `Status: proposed` proposals, then `handoff.md` last and only if under 7 days old. The board wins on any disagreement.
- **Promotion and expiry:** unchanged path (observation, cheap reproduction, `proposed` board task, human `ready`, narrow edit to the owning skill with a check). Observations older than 30 days are promoted or deleted. The whole folder is read at each `VERSION` bump.
- **Folder name:** the review suggested the on-disk name `lab-files` (spaces are a quoting hazard). The maintainer first chose "Sensei Ewok Lab Files", then renamed it to `cf-lab-files` the same day, together with the repos (`hq` to `cf-lab`, `research` to `cf-research`, `skills` to `cf-skills`). The name no longer has spaces; the setup test still runs from a root folder that contains a space.
- **Still open:** the deterministic check for the folder (board T-0056); whether hq's deny rules for `Read(.env)` apply to files outside the hq root in a multi-root workspace (unverified); backup, since the folder has no history.
