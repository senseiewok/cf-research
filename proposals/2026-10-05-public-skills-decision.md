# Using Anthropic's public skills: the council's decision

Status: proposed, written 2026-10-05. A council of three (Fable, the local Qwen3.8 27B in thinking mode, and the controlling agent, each answering the same briefing without seeing the others) was asked how the lab's public repos should use `github.com/anthropics/skills`. Their rankings, best first:

| Member | Ranking |
| --- | --- |
| Controlling agent | keep the clone private, point to skills by name, install per machine, make our own skills installable, nothing vendored |
| Fable | keep the clone private, install on a maintainer machine, point to skills by name, make our own skills installable |
| Qwen | point to skills by name, install on a maintainer machine, make our own skills installable, keep the clone private |

All three ranked vendoring copies into a public repo and a git submodule last. They differed on one point: Qwen would add a one-line pointer to the public repos now, while Fable and I would add none until a trial shows a skill earns it. Qwen's own risk list argued against it (a pointer can quietly become policy), so the decision follows Fable and the controlling agent.

## Decision

1. **Now** (the part about the public repos is superseded by the update below): the clone stays a read-only folder in the maintainer's private admin workspace (commit 8a1541c, 2026-09-28). Nothing in the public repos changes.
2. **Trial first, small and measured** (board T-0075): `frontend-design` on one private page with and without it, judged by checks that run (accessibility, no new dependency, verified copy untouched); `skill-creator`'s with-and-without evaluation idea on one existing lab skill that has a deterministic check, counting only gains the check confirms. Read every script before running any of it, and run nothing unattended. Measure what the skill descriptions cost in context.
3. **After a positive trial only** (superseded by the update below, which allows the pointer now): one labelled-optional pointer in the lab skill that benefits (name, link, pinned commit, "third party, not reviewed by the lab"), and a short evaluation step written in the lab's own words with attribution. Never copied text, never a vendored skill, never a submodule.
4. **Our own skills as a plugin** (board T-0076): add `.claude-plugin/marketplace.json` to the skills repo so others can install the lab's skills with one command. It serves Claude Code users only; GitHub Copilot still reads `.claude/skills`.

## Update, 2026-10-05, later: the maintainer's decision

After the council, the maintainer decided that the lab's skills are layered on top of Anthropic's in the lab's own repos. This takes the step the council had held back until after the trial (the pointer, option O5) and keeps its other limits: nothing is copied or vendored, the base is named with its source, licence and reviewed commit and labelled third party, and each lab skill must work alone for someone who has not installed the base. The pairings are listed in `THIRD_PARTY_SKILLS.md` in the control repo and the rule is section 11 of `ai-provider-compatible-skills`. The trial (board T-0075) still decides whether to adapt any of the base skills' ideas into the lab's own wording.

## Facts the council corrected or left open

- The clone does have a `.claude-plugin/marketplace.json`; the controlling agent's earlier claim that it had none was wrong (the listing hid dot-folders). It is a marketplace of bundles: `example-skills` holds twelve skills at once (algorithmic-art, brand-guidelines, canvas-design, doc-coauthoring, frontend-design, internal-comms, mcp-builder, skill-creator, slack-gif-creator, theme-factory, web-artifacts-builder, webapp-testing), `document-skills` holds the four proprietary document skills, and three more are single-skill plugins. So a plugin install cannot be selective: choosing `skill-creator` means loading all twelve descriptions every turn, including one with no licence file. That is the strongest argument for trying things by reading the clone, and for adapting an idea into a lab skill instead of installing a bundle.
- Licences, verified on the files: Apache-2.0 for the skills of interest, Anthropic's proprietary licence for docx, pdf, pptx and xlsx, no licence file for doc-coauthoring.
- Whether Claude Code discovers skills in an added workspace folder is unsettled: the documentation summary says no, while this session saw skills under an added folder's `.claude/skills` appear at once. One reproducible check is part of T-0075; do not build on either answer yet.
- The research agent read most skills through a summariser and did not run any; the fit ratings are its opinion until the trial.

## Update, 2026-10-05, evening: how the skills repository was built, and what the open questions became

- **Skills in an added folder do load.** In the maintainer's Claude Code session `cf-skills` was an additional working directory, and the skills under its `.claude/skills/` (`cf-evidence-loop`, `site-seo-review`) appeared in the session's list of available skills as soon as they existed. This settles the open question above for Claude Code in this setup; it was not repeated in other agents or versions.
- **Layout.** At the maintainer's request the public skills repository keeps its skills in `.claude/skills/` (the folder Claude Code, Copilot and Cursor read from a project), so a clone works as a skills folder at once. Other material lives in ordinary folders: `docs/` (per-agent install notes and copy-paste prompts, each path dated and labelled documented or partial), `scripts/` (`install_skill.py`, `make_index.py`, `check_repo.py`), `template/`, `spec/`, and a generated `skills-index.json` that lists what each skill can do (network, runs code, needs Python, optionally a browser) so a person is told before they install. The plugin marketplace is optional and lists each skill as its own plugin.
- **Nothing of Anthropic's is vendored or copied.** `docs/recommended.md` lists four of their skills that the lab tried, with licences and read-first notes, and the single-skill `gh skill install` command that GitHub documents; Anthropic's own plugin bundles load about a dozen skills at once, which the single-skill route avoids. The four local copies in the control repository stay git-ignored trial installs (board T-0075).
- **Not yet verified.** Nothing was run in a real client (board T-0082); `gh` is not installed on the maintainer's machine, so the `gh skill install` commands are built from GitHub's documentation and have not been run.
