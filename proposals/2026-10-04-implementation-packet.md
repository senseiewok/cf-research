# Task packet: apply the 2026-10-04 registry-access research to the workspace

Status: superseded by `2026-10-04-finalize-packet.md`, which includes it. Kept for the record.

Paste everything below the line into the orchestrating agent in `hq`. It is written to the `ai-loop-council` packet contract: goal, grounding files, constraints, bounded steps, deterministic acceptance checks. The orchestrator owns the diff; the local worker may draft individual steps.

---

## Goal

Apply the reviewed changes in `../research/proposals/2026-10-04-hq-skill-updates.md` to the `hq` skills, wire the new `research/sources/catalog.yaml` into the workspace's agent instructions as the single record of source-access permission, and leave both repos in a committable state. Do not extend scope beyond what the proposal file names.

## Read first, in this order

1. `hq/AGENTS.md` (workspace rules; the "Improve the existing skills" bullet governs this task)
2. `hq/.claude/skills/security-browsing/SKILL.md` and `hq/.claude/skills/cf-research-context/SKILL.md` (edit targets)
3. `../research/proposals/2026-10-04-hq-skill-updates.md` (the exact text to apply, with evidence)
4. `../research/sources/README.md` and `../research/sources/catalog.yaml` (what the skills will point at)
5. `../research/landscape/registry-data-access.md` (the evidence behind every claim; do not re-derive it)

## Constraints

- **Apply, do not author.** Use the proposal's text verbatim. If a passage seems wrong, stop and report it; do not improve it in place.
- **No new facts.** Every claim in the inserted text already carries a source and a check date. Do not add numbers, URLs, or dates from model memory.
- **No fetching.** This task needs no web access. The one permitted exception is the verification step below, which runs a local script that fetches from `pr.ecfs.eu` only; the script refuses anything the catalog does not permit.
- **Never touch `www.cff.org`, `databases.lovd.nl`, or `cftr2.org` programmatically.** Not with `requests`, not with Playwright, not with any client. The reasons are in the evidence file; they are not up for re-evaluation here.
- **Protected paths.** `hq/.claude/skills/*/SKILL.md` edits require a human to approve the diff before commit (Elevated review tier: agent instructions). `hq/.claude/settings.json` is out of scope.
- **Frontmatter stays valid.** `name` must match the folder; keep `license` and `compatibility` unchanged.
- **No PHI, no secrets, no local machine details** in any edited file.

## Steps

1. **`security-browsing`** — insert the tested-results table and the catalog pointer after section 4, rule 5, exactly as proposal item 1 specifies. Do not delete the existing rule text.
2. **`cf-research-context`** — (a) replace the intro sentence of "Data and community sources" and add the ECFSPR bullet per proposal item 2; (b) insert the "Registry numbers" subsection after "Disease basics" per proposal item 3, preserving its observed/inferred labels.
3. **`ai-loop-council`** — add the "Variant: value extraction from a document" block under the task packet template, per proposal item 4.
4. **`model-onboarding`** — add proposal item 5 as a clearly labelled *design only* paragraph in section 4. Build no fixture.
5. **Catalog smoke test** — from `research/`, run:
   ```
   pip install -r tools/sources/requirements.txt
   python tools/sources/fetch_sources.py --manual
   python tools/sources/fetch_sources.py --only ecfspr-adr-2022
   python tools/sources/fetch_sources.py --verify
   ```
6. **`.github/loop-orchestrator/loop.py`** — add a hard failure at the top (`raise SystemExit("not implemented: this prototype reports success unconditionally")`) so it cannot be mistaken for a working verifier. Leave the file otherwise untouched; update its mention in `AGENTS.md` to say it now exits non-zero.
7. **Commit messages** — one commit per repo. `research`: "Add source catalog, permitted-source fetcher, and HQ skill proposals". `hq`: "Apply registry-access research to browsing, CF context, loop, and onboarding skills". Do not push; the human pushes.

## Acceptance checks (verifier decides, not a model)

- [ ] `python -c "import yaml; yaml.safe_load(open('research/sources/catalog.yaml'))"` exits 0.
- [ ] Every `SKILL.md` edited still begins with `---`, parses as YAML frontmatter, and its `name` equals its folder name. A short PowerShell or Python check that asserts this for all four files is acceptable evidence; a model saying "frontmatter looks fine" is not.
- [ ] `fetch_sources.py --manual` output lists only `access: manual` catalog entries (the Cystic Fibrosis Foundation reports and CFTR2) and no `pr.ecfs.eu` URL.
- [ ] `fetch_sources.py --only ecfspr-adr-2022` downloads one file; `--verify` then prints `1 match, 0 changed, 0 missing` and exits 0.
- [ ] `git -C research status --porcelain` shows no file under `sources/downloads/` and no `*.part`.
- [ ] `grep -ri "playwright" hq/.claude/skills/security-browsing/SKILL.md` matches only the sentence forbidding it for CFF, plus any pre-existing mentions.
- [ ] `python hq/.github/loop-orchestrator/loop.py` exits non-zero.
- [ ] `git diff` in `hq` touches only: the four `SKILL.md` files, `AGENTS.md`, `.github/loop-orchestrator/loop.py`. Anything else is a scope violation; revert it.

## Escalate to the human if

- Any proposal passage conflicts with the current skill text in a way the proposal did not anticipate (the skill changed since 2026-10-04).
- A check above fails twice after one bounded repair attempt.
- Any step would require network access beyond `pr.ecfs.eu` or a package not in `tools/sources/requirements.txt`.
- The worker proposes to "clean up" anything outside the listed files.

## Report back with

Command output for every acceptance check, the two `git diff --stat` summaries, and a one-line list of anything skipped and why. No prose summary of what the skills now say; the human will read the diff.
