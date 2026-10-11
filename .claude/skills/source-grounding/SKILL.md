---
name: source-grounding
description: How to ground every public sentence that states a fact, a definition, a term or a number before it is published, covering the grounding ladder (computed from data, a registry statement checked verbatim, a definition quoted from an open-access primary source, a lab position labelled as such), the answer tiers T0 to T3, what each checker proves and does not prove, the licence rule for quotes, the guard against invented claims, how a local model may help safely, and Spanish (es-MX) terms. Use before writing or reviewing any page text, caption, gloss, definition, claims file or note that others will read, and before asking a model to find or extract a definition.
license: CC0-1.0
compatibility: Works with any Agent-Skills-spec-compatible tool that reads .claude/skills (Claude Code, GitHub Copilot; see AGENTS.md in cf-lab). The checkers it names need Python 3.9 or newer; the model steps need PowerShell 7 and a local Ollama model. A local worker model does not load skills; it gets pasted packet text from references/local-worker-packets.md.
---

# Source grounding

Every public sentence that states a fact, a definition, a term or a number rests on something a stranger could check. This skill is the catch-all: the ladder, the checkers and the guards. Specifics seen in practice are in `references/lessons.md` until promoted. The disease rules and the tiers are owned by `cf-research-context` in cf-lab; this skill applies them.

## How each agent gets this skill

| Agent | How |
| --- | --- |
| Claude Code | Start in `cf-lab` and add this repository with `--add-dir ../cf-research`, or `/add-dir`; AGENTS.md in cf-lab says skills from an added folder reload live. `permissions.additionalDirectories` alone gives files, not skills |
| GitHub Copilot | AGENTS.md in cf-lab says Copilot reads `.claude/skills/`. Whether it reads that folder in every folder of the multi-root workspace was not tested here: see AGENTS.md |
| A local worker model | Loads no skill. Paste a packet from `references/local-worker-packets.md` |

## The ladder

Use the highest rung that fits. A sentence that fits none is not published.

| Rung | What the sentence rests on | How it is checked |
| --- | --- | --- |
| 1. Computed | A count or value computed from data by a script | The script recomputes it; drift fails (for the atlas, `check_atlas.py`) |
| 2. Registry statement | A sentence the registry says about itself or its data | Verbatim quote in a claims file (`tools/trial_atlas/registry_claims.json`), checked with `tools/claims/check_claims.py` through `check_registry_claims.py` against the saved page text |
| 3. Definition | A definition quoted from a primary open-access source | The grounding tool finds it and verifies the quote as an exact substring of the article (in progress: see below) |
| 4. Lab position | The lab's own reading or choice | Labelled as the lab's position or "our inference" in the sentence itself |

## Tiers

As in `cf-research-context`: T0 a summary or memory, unverified; T1 a quote from a primary source with its id and date; T2 T1 plus a script check and a different model's check, each flag confirmed on the source; T3 T2 plus a person's review. Public claims need T2; scientific copy on a site and anything that could affect care need T3. Say the tier.

## What each checker proves, and does not

| Checker | Proves | Does not prove |
| --- | --- | --- |
| `tools/claims/check_claims.py` | Each quote is exact; each number is in the quote or computed from it; widening words and absence claims carry a scope; inferences carry a basis | That the quote supports the sentence |
| `tools/trial_atlas/check_registry_claims.py` | Every registry sentence the tools emit is covered by claims, and (with saved pages) every quote is on the registry's page | That the claim is a fair reading of the page |
| `check_atlas.py` model-tag check | A model's quote is an exact substring of the entry | That the class is right |
| `tools/sources/check_source_overlap.py` | No quotation over 8 words, no long copied runs | That a shorter quotation is fair or true |
| `tools/sources/check_numbers.py` | Every number in a draft exists in the evidence | That it sits on the right claim |
| The grounding tool | The quoted words exist in an open-access article | That the definition is right |

A verified quote proves only that the words exist. A different model or a person still reads each quote against its sentence.

## The grounding tool (in progress)

A tool for rung 3 is being built for the atlas (T-0130); its README will be the single source for how to run it. The method: searches for each term go to an open-access literature index with an open-licence filter; a model proposes a short definition quote; a script keeps only a quote that is an exact substring of the article and holds the term; the licence of every hit is checked again. Its search fields are recorded as "to verify on the first real run" (see `references/lessons.md`).

## Licence rule for quotes

- Quote only from openly licensed text whose licence the tool or a person has confirmed.
- Article text and model outputs stay in private folders, never in a repository.
- Public text carries short quotes only, within the lab's quotation limit, with the source linked.

## The guard against invention

1. A model proposes; a script verifies (exact quote, numbers); a different model or a person reads each quote against its sentence.
2. Widening words ("only", "all", "same", "no longer", "first") need a written scope of what was searched.
3. A claim that something is absent needs a scope record.
4. Where the source is silent, write "not stated". Never fill a gap from memory.

## Using a local model safely

- Wrap every untrusted text in `<untrusted_page>` tags with the data-boundary wording (security-browsing, section 1). Give the model no tools; ask for JSON only, against a schema.
- Verify every quote as an exact substring by script before reading the reply as prose.
- Name the model profile by file (for the lab's worker, the 64K profiles in `model-qwen3-8-27b`). With a profile set, do not rely on `LOCAL_WORKER_MODEL`: it replaces the profile's model name. The context size comes from the profile, not from a model alias name.
- End every packet that asks for facts with the closing lines in `ai-loop-council`.

## Spanish terms

Spanish means Mexican Spanish (es-MX). A Spanish term needs a quote that shows its use in that article. The maintainer reviews every Spanish term before it is published. Registry quotes and instrument names stay as written.

## Lessons: how specifics become rules

Specifics live in `references/lessons.md` as dated observations (date, what was seen, tier, a "reproduce with:" line, status: observed, proposed-rule or promoted). An observation becomes a rule only through a `proposed` board task, a person setting it `ready`, and a narrow edit to the owning skill with a check (AGENTS.md, "Local files and memory" and "AI loop and delegation"). Delete the observation once promoted. Review or remove an unpromoted observation older than 30 days. Never promote a model's opinion into policy.

## Verification

- `python -m unittest -v` in `tools/claims` and in `tools/trial_atlas` (no network).
- Before publishing: the claims file passes `check_claims.py`, the draft passes `check_source_overlap.py --quotes-only`, and a different model or a person has read each quote against its sentence.
