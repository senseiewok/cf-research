# Task packet: publish `cf-evidence-loop` and the supporting research changes

Paste everything below the line into the orchestrating agent in `hq`. It follows the `ai-loop-council` packet contract. It supersedes `2026-10-04-finalize-packet.md` for the `skills` and `research` repos; the `hq` skill-edit steps in that packet still apply and are referenced, not repeated.

---

## Goal

Bring `skills/cf-evidence-loop` through the `skills` repo's admission requirements and to a reviewable pull request; bring the `research` repo (source catalog, fetcher, proposals, pointer) to a PR; apply the `hq` skill edits already specified. Publication means a human merges; you prepare.

## Read first, in this order

1. `hq/AGENTS.md` and `hq/tasks/README.md`; then `hq/tasks/BOARD.md`. You are executing **T-0003, T-0004, T-0022, T-0037, T-0038, T-0048** and preparing PRs. Set them `in_progress`; a human sets `done`.
2. `skills/README.md`: the five admission requirements. They are the acceptance test for the skill.
3. `skills/cf-evidence-loop/SKILL.md`, then `references/ARCHITECTURE.md`, then `catalog.yaml`.
4. `hq/.claude/skills/ai-provider-compatible-skills/SKILL.md` and `security-git/SKILL.md`.
5. `research/proposals/2026-10-04-hq-skill-updates.md` and `2026-10-04-implementation-packet.md` for the `hq` edits.
6. `research/proposals/2026-10-04-mcp-and-api-survey.md`: context only; its tasks are not in this packet.

## Constraints

- **The license question is open (T-0001, P0).** The skill's frontmatter says `CC0-1.0` to match existing skills, but this is the first package with code. Do not pick a license. Stop at step A4 and ask the human: CC0-1.0 for everything, or MIT/Apache-2.0 for code with CC0/CC-BY for prose. Record the answer on the board and apply it to all four repos.
- **Apply, do not author.** Skill text for `hq` comes verbatim from the proposal. The `cf-evidence-loop` code is reviewed by you, not rewritten: fix a defect you can demonstrate with a failing test; do not restyle.
- **No new providers** in this packet. The one new source, Data Commons (section D), enters the catalog as `unverified` through its admission review; its provider is a later task.
- **Network conduct:** `skills/cf-evidence-loop/references/NETWORK-RULES.md` binds every run in this packet. Run each live check once. A `429`, `403`, or cooldown line is a result to report, never a reason to retry, change clients, or wait-and-loop. Space the live checks at least a minute apart.
- **Network calls permitted:** only those the skill's `catalog.yaml` permits, and only through the skill's own client, for the live checks named below. Never any client against `www.cff.org`, `databases.lovd.nl`, or `cftr2.org`.
- **Never commit** anything under `sources/downloads/`, any `.jsonl` ledger, `__pycache__`, or `.pytest_cache`.
- **No PHI, secrets, machine details, local user paths, or session artifacts** in any file or PR body. `grep -rn "C:\\Users" skills research` must return nothing in tracked files.
- **Do not push, do not merge, do not move board tasks to `ready` or `done`.**

## Steps

### A. `skills` repo

1. Branch `feat/cf-evidence-loop-v0.1`.
2. Run the offline suite: `python -m pytest -q tests` from `skills/cf-evidence-loop`. Expect 37 passed, including tests/test_conduct.py (one per network rule).
3. Run the live checks from a directory outside the repo, with no `PYTHONPATH` and no `EVIDENCE_CATALOG` set:
   - `python <path>/scripts/evidence_cli.py retractions 10.1016/S0140-6736(97)11096-0` returns exactly two notices, types `correction` (2004) and `retraction` (2010).
   - `python <path>/scripts/evidence_cli.py approval TRIKAFTA` returns `NDA212273`, `original_approval=2019-10-21`.
   - `python <path>/scripts/evidence_cli.py sources` lists nine permitted sources and none of `lovd-cftr`, `cffpr-adr-2022`, `cffpr-technical-supplement`.
4. **Stop and ask the human the license question** (see constraints). Apply the answer: `LICENSE` file at the `skills` repo root, matching `license:` in the frontmatter.
5. Admission checklist against `skills/README.md`, with evidence for each line:
   - Useful outside the workspace: step 3 ran from outside the repo with the skill's own catalog.
   - Canonical `SKILL.md` with matching `name`: a script asserts frontmatter parses, `name == cf-evidence-loop`, description 1-1024 chars, compatibility <= 500 chars.
   - License and provenance declared: step 4; `ARCHITECTURE.md` names the upstream data sources; `catalog.yaml` carries check dates.
   - Runnable checks: step 2 output.
   - One maintained source: `research/tools/evidence/README.md` is a pointer only; `grep -rn "def parse_updates" research skills` finds exactly one definition.
6. Register the package in `skills/README.md` (already present; change its status line from "awaiting admission review" to the outcome).
7. Write the PR description to `skills/PR-2026-10-04.md` using the template at the end.

### B. `research` repo

1. Branch `feat/source-catalog-evidence-pointer`.
2. `python -c "import yaml; yaml.safe_load(open('sources/catalog.yaml'))"` exits 0.
3. Fetcher smoke test: `python tools/sources/fetch_sources.py --manual` lists only Cystic Fibrosis Foundation entries; `--only ecfspr-adr-2022` then `--verify` prints `1 match, 0 changed, 0 missing`. If `tools/sources/test_fetch_sources.py` exists, run it and report.
4. Confirm `tools/evidence/` contains only `README.md` (the pointer).
5. `git status --porcelain` shows nothing under `sources/downloads/`, no `*.part`, no `*.jsonl`.
6. Write the PR description to `research/PR-2026-10-04.md`.

### C. `hq` repo

1. First: `git status` will show roughly 20 modified files that predate this work. **Stop and ask the human** whether to commit, stash, or discard them. Exclude them from this PR.
2. Branch `chore/registry-access-skills-and-board`.
3. Apply the four skill edits per `2026-10-04-implementation-packet.md` steps 1-4, and the `loop.py` hard failure per its step 6.
4. In `AGENTS.md`, under "Skills / Domain context", add one line pointing at `../skills/cf-evidence-loop` as the evidence-record tool agents use before citing a number, approval, trial, or paper. Do not copy the skill into `hq/.claude/skills/`; your rule is one maintained source and an install/link mechanism decided before moving. Record the link mechanism decision as a new board task if none exists.
5. Add `tasks/render-board.ps1 -Check` to the repo's existing validation step; if none exists, note it under "Validation expectations" in `.github/copilot-instructions.md`.
6. Write the PR description to `hq/.github/PR-2026-10-04.md`.

### D. Data Commons MCP admission review (board T-0048)

Purpose: decide whether Google's Data Commons MCP server may become the lab's population-denominator source. This is a review, not an integration. Nothing here is added to `cf-evidence-loop` until the review passes and a human says so.

Facts established 2026-10-04 from the PyPI record, to verify against the repository: package `datacommons-mcp` 1.4.0, Apache-2.0, Python >=3.11 <3.14, source `github.com/datacommonsorg/agent-toolkit`, requires `DC_API_KEY` (free, apikeys.datacommons.org), serves over stdio and streamable HTTP, depends on `fastmcp==3.4.2` and `google-cloud-storage`. Data Commons hosts were unreachable from the sandbox where this packet was written (allowlist, not Google), so every live step below runs on the lab machine.

1. **Provenance** (`model-onboarding` section 1 applied to a server): confirm the GitHub organisation is `datacommonsorg`; pin release 1.4.0 by exact version; read `pyproject.toml` and any install or post-install scripts; list every runtime dependency and note that `google-cloud-storage` is present and ask what it is used for. Install into a fresh virtual environment, never the system Python.
2. **Terms and rate limits**: read the Data Commons API terms and the API-key documentation. Record `terms_url`, any stated rate limit, and the data licence of the underlying statistics (Data Commons aggregates many sources with their own licences; the provenance of each statistic is what the lab needs, not a blanket statement). Add a `datacommons` entry to `research/sources/catalog.yaml` with `access: api`, `auth: api_key_required`, `max_rps: 1.0`, `terms_url`, and `claim_label: unverified` until step 4 passes.
3. **Key handling**: the key lives in the lab's local, git-ignored `.env` or OS credential store (`security-git`). It never appears in a task packet, a board note, a PR, or a worker prompt. `grep -rn DC_API_KEY hq research skills` must match only documentation sentences, never a value.
4. **Runtime review** (`security-runtime` MCP checklist): run the server over **stdio only**, not HTTP, attached to one agent client with fetch approvals on. Issue exactly three queries, by hand, a minute apart: (a) "population of the United States, latest", (b) "population of Washington state, latest", (c) "cystic fibrosis" as an indicator search, then (d) "sickle cell disease" as an indicator search. Record each response's source and date fields verbatim.
5. **Decision**: write `research/proposals/<date>-datacommons-admission.md` with: pass/fail per step, the four responses, and a one-paragraph recommendation. The expected finding is that (a) and (b) return population figures with a US Census source and date, and that (c) and (d) return nothing CF- or SCD-specific; if (c) or (d) does return a disease indicator, record its source precisely, because that would change the registry project's scope.
6. **If it passes**: append a board task for a `datacommons` provider in `cf-evidence-loop` that follows NETWORK-RULES.md and returns Evidence records whose `limitations` name the underlying source and vintage. Do not write the provider in this packet.

Escalate and stop if: the install runs a script you did not read; the key cannot be kept out of the repos; a query returns anything that looks like individual-level data; the terms restrict redistribution of the statistics in a way that affects the registry table.

## Acceptance checks (verifier decides; paste real output)

- [ ] `pytest`: `37 passed` in `skills/cf-evidence-loop`.
- [ ] Live retraction check output showing two notices with the stated types and years, and the stderr accounting line (`requests: N | by host: ... | cooldown: none`). If it reports `rate_limited` or a host in cooldown, record that as the result and stop; do not re-run within the hour.
- [ ] `EVIDENCE_DRY_RUN=1` run of any command prints planned URLs, `requests: N`, and opens no socket (verify with the request log empty of real statuses).
- [ ] Live approval check output showing `NDA212273` and `2019-10-21`.
- [ ] `sources` output: nine lines, none of the three refused ids.
- [ ] Frontmatter assertion script output: pass.
- [ ] `grep -rn "def parse_updates" research skills`: exactly one hit, under `skills/`.
- [ ] `LICENSE` present at the root of `skills`, `research`, `hq`, `.ai`, each naming the license the human chose.
- [ ] `yaml.safe_load` on `research/sources/catalog.yaml`: exit 0.
- [ ] Fetcher `--manual` and `--verify` output as specified.
- [ ] `pwsh -File hq/tasks/render-board.ps1 -Check`: `BOARD.md is current`.
- [ ] `python hq/.github/loop-orchestrator/loop.py`: non-zero exit.
- [ ] `git -C skills status --porcelain` and `git -C research status --porcelain`: no `__pycache__`, `.pytest_cache`, `*.jsonl`, `sources/downloads/`, `*.part`.
- [ ] `grep -rn "C:\\Users" skills research hq/tasks hq/AGENTS.md`: no hits.
- [ ] `git diff --stat` for each repo touches only the files this packet names.
- [ ] Data Commons: `pip show datacommons-mcp` shows version 1.4.0 in a fresh virtual environment; the catalog has a `datacommons` entry with `terms_url` and `max_rps`; `grep -rn DC_API_KEY hq research skills` matches only prose; the admission note exists with four recorded responses, each carrying a source and a date; `git status` in every repo shows no `.env`.

## Escalate to the human, and stop, if

- The license question (always; it is a required stop).
- The `hq` pre-existing modifications (always; required stop).
- Any live check returns `blocked` or `rate_limited`. Report it as the publisher's answer; do not retry with a different client.
- A proposal passage conflicts with current skill text.
- A test fails twice after one bounded repair.
- You want to add a provider, source, or feature. Append a board task instead.

## Report back with

1. The three `git diff --stat` outputs.
2. Real output for every acceptance check, in order.
3. Paths of the three PR description files.
4. One line per board task touched: id, status, the acceptance output that justifies it.
5. One line per thing skipped or appended to the board, with the reason.

No prose summary of what the skill does. The human reads `SKILL.md`.

---

## PR description template

```
## What
<one paragraph; board task ids this PR completes>

## Why
<two or three sentences; link the proposal files applied>

## Evidence
<each acceptance check and the deciding line of its real output>

## Human review items
<things a model must not decide: license choice, skill wording, scope triage>

## Provenance
Proposed by: <created_by from the board tasks>
Applied by: <your model_id>, <date>
Review tier: <routine | elevated | full>

## Not in this PR
<board task ids deliberately left proposed>
```
