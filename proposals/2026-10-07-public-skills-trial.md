# Anthropic's public skills: the trial (board T-0075)

Status: proposed, written 2026-10-07. The maintainer decides, per skill, whether to adopt, adapt or drop. Nothing here changes a lab skill. It follows the plan in `2026-10-05-public-skills-decision.md`: try `frontend-design` with and without, test `skill-creator`'s with-and-without idea on a lab skill that has a deterministic check, read every script before running it, and measure what the descriptions cost.

Source: github.com/anthropics/skills, commit 8a1541c (2026-09-28), Apache-2.0 for the skills below. They are third-party instructions, read before use and run nothing unattended.

## What was run

| Trial | Setup | Judged by |
| --- | --- | --- |
| `frontend-design`, with and without | The local worker (Qwen3.8 27B, fast mode) wrote one HTML page from the same task, five seeds per arm. The "with" arm had the skill body prepended | `tools/skill_trial/check_page_basics.py`: seven checks (dependencies, lang and landmarks, contrast, mobile at 360 px, motion, verified copy word for word, console) |
| `ascii-art` (a lab skill), with and without | The same loop: a lantern for a README, five seeds per arm. This tests the with-and-without idea from `skill-creator`; it is not a test of an Anthropic skill | `tools/skill_trial/check_ascii_art.py`: wraps `check-ascii.py --strict --max-width 36` and requires at least one art block, because that checker passes a file with no art at all |
| `skill-creator` | Read in full. Of its scripts only `quick_validate.py` was run, on lab skills. Its `claude -p`, server and port-killing scripts were read and not run | Reading, and the one script's output |
| Discovery | `tools/skill_trial/check_skill_discovery.py` | Claude Code 2.1.292, Windows, two runs per arm |

The checker was tested before it judged anything: `test_check_page_basics.py` holds a good page, a good page with a guarded animation, a good page with an inline SVG, and 17 bad variants, each the good page with one change, and prints `VERIFIED` only when every bad variant fails the check it should. A mistake found along the way (a `url(#...)` inside a `data:` URI was counted as an outside load) was fixed, tested, and the saved pages were re-judged without a new model call.

## `frontend-design`

| Arm | Mean checks passed (of 7) | Pages passing all 7 | Prompt tokens | Mean output tokens | Mean seconds |
| --- | --- | --- | --- | --- | --- |
| Without | 5.6 | 1 of 5 | 242 | 2,928 | 20.8 |
| With | 6.4 | 3 of 5 | 2,178 | 2,554 | 18.8 |

What failed, per page:

- Without: seed 12 changed one of the three verified sentences and had an animation that runs forever; seed 13 had no `main` landmark and low contrast (3.22); seed 14 had 9 text elements over a gradient or image that could not be measured; seed 15 had contrast 4.34 and two endless animations.
- With: seed 11 had 2 text elements over a gradient or image that could not be measured (the lowest measured ratio, 5.16, passes); seed 13 had no `main` landmark and low contrast (3.13).

Read this as a direction, not a result. Five pages per arm, one task, one model, and a checker written for this trial. The checks measure mechanical quality (accessibility, no outside loads, copy untouched), not what the skill is for, a distinctive look. The "with" arm was no slower and wrote fewer tokens, but it carries 1,936 more prompt tokens on every call. A reviewer has not yet looked at the ten pages for taste; only the checks and the controlling agent's own screenshots were used.

What in the skill looks worth keeping, from reading it: spend boldness in one place, cut decoration that serves nothing, use motion sparingly and only for a reason, do not number or label content that is not a sequence, a quality floor (mobile, focus, reduced motion, contrast) built in without announcing it, and a two-pass plan-then-critique process. What does not fit the lab: it pushes toward an opinionated, risk-taking visual identity, which the calm, sober audience brief for this lab works against. The lab voice and verified-copy rules win where they differ.

## `ascii-art` (lab skill)

| Arm | Pages passing both checks | Prompt tokens |
| --- | --- | --- |
| Without | 1 of 5 | 139 |
| With | 4 of 5 | 5,473 |

Without the skill, seeds 13, 14 and 15 had art with no description after the block (the accessibility rule), and seed 12 had no fenced block that the checker could find. With the skill, seed 13 returned unfenced Braille dot art with a description after it, so the checker found no block. The skill earns its cost on this task, mostly by teaching the fenced block and the description-after-the-block convention, but five runs on one task only show that the check can tell the arms apart.

## Cost in context

Measured with the local worker's tokenizer (Qwen), so other models count differently. The description is loaded every turn when the skill is available; the body is loaded when the skill is used.

| Skill | Description (every turn) | Body (when used) |
| --- | --- | --- |
| `frontend-design` | 41 | 1,916 |
| `skill-creator` | 65 | 7,466 |
| `ascii-art` (lab) | 162 | 5,312 |

Across the 23 skills measured (lab and third party), descriptions range from 41 to 162 tokens per turn. Anthropic's own plugin bundle `example-skills` loads twelve skills at once and cannot be installed selectively (see the decision memo), which is why the trial used hand-copied, git-ignored folders.

## `skill-creator`

Read: its SKILL.md (a method for writing and testing skills) and its scripts. Findings from reading, not from running, except where stated.

- Its with-and-without evaluation is the idea the `ascii-art` trial above reproduces in the lab's own terms: same task, with and without the skill, judged by a check that runs and counted only where the check confirms. The trial adds the vacuous-pass guard (a check must find something to check) and a test that the check can fail.
- Its scripts call `claude -p` on the user's own login, write a temporary command file into the project and serve a review page on 127.0.0.1. Not run, and not needed.
- `quick_validate.py` was run on lab skills. On Windows it failed with an encoding error until `PYTHONUTF8=1` was set, because it reads the file without naming an encoding. From reading the code, its frontmatter pattern also requires `\n` line endings, so a CRLF file would not match; this was not run. Reported as an observation, not fixed upstream.
- Three of its rules are not in the public skills repo's own check (`../cf-skills/scripts/check_repo.py`, which already checks the name pattern, the 64-character name limit and the 1024-character description limit): no frontmatter keys outside `name`, `description`, `license`, `allowed-tools`, `metadata` and `compatibility`; no angle brackets in the description; `compatibility` at most 500 characters.

## Discovery

`check_skill_discovery.py`, run in this session on Claude Code 2.1.292 (Windows), cheapest model, only the `Skill` tool enabled: without `--add-dir` the probe skill was not listed in either run; with `--add-dir` pointing at a folder holding it, it was listed in both runs. With every tool disabled the model reports no skills list at all, so a control run in that state would prove nothing; that is why the check enables `Skill`. This settles the open question in the decision memo for this version and setup only.

## Recommendation (for the maintainer)

| Skill | Suggested decision | Why |
| --- | --- | --- |
| `frontend-design` | Adapt, do not install. Add restraint rules (one bold element, motion only with a reason, a quality floor) to the lab's own web guidance in the lab's words, with attribution | Positive direction on mechanical checks at small n; the skill's main push (a distinctive, risk-taking identity) does not fit the audience; 1,936 prompt tokens per use |
| `skill-creator` | Adapt, do not install. Add the three missing frontmatter rules to `check_repo.py` with a test for each, and keep the with-and-without method as described above | The with-and-without idea is reproduced above in the lab's terms; its scripts do model calls, serve a page and kill processes; one script fails on Windows without a fix; 7,466 tokens to load |
| `algorithmic-art`, `theme-factory` | Not trialed. Leave as recorded in `THIRD_PARTY_SKILLS.md` | No evidence either way |

If adapted, each change is its own small pull request with its own test, and each factual or design claim is in the lab's words. Adoption of a lab skill from this trial stays a human decision.

## Not shown

- Any gain on a task other than a single page and a single lantern, or with another model, or at a larger n.
- That the pages are better to look at; the checks do not measure taste.
- Behaviour in Copilot or another agent; only Claude Code 2.1.292 was checked.

Reproduce with: `python tools/skill_trial/test_check_page_basics.py` (expects `VERIFIED`); `python tools/skill_trial/run_arms.py --help` for the with-and-without runner; `python tools/skill_trial/check_skill_discovery.py` (two small model calls on your own login).
