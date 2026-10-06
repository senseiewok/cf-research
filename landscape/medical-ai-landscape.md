# Medical AI and agent-skills landscape

External projects and references relevant to building CF-focused AI agents, skills, and MCP servers. Every entry is labelled per the [README](../README.md#how-to-read-claims) claim legend.

**Status key:** `[verified]` confirmed live · `[unverified]` plausible, not confirmed · access date shown for verified items.

> Snapshot date: **2026-10-01** (star counts and activity change; re-verify before citing).

## Agent skill repositories

| Repo | What it is | Status |
| --- | --- | --- |
| [`FreedomIntelligence/OpenClaw-Medical-Skills`](https://github.com/FreedomIntelligence/OpenClaw-Medical-Skills) | "Largest open-source medical AI skills library for OpenClaw." Python. **~3.0k★**, active as of Jul 2026. | `[verified]` (GitHub search, 2026-10-01) |
| [`K-Dense-AI/scientific-agent-skills`](https://github.com/K-Dense-AI/scientific-agent-skills) | Scientific agent skills; part of the `K-Dense-AI` org (also `mimeographs`, ~129★). Widely forked/packaged by third parties. | `[verified]` org + forks (2026-10-01) |
| [`xjtulyc/MedgeClaw`](https://github.com/xjtulyc/MedgeClaw) | An open-source AI research assistant for biomedicine, by its own description (RNA-seq, drug discovery, clinical analysis). Built on Claude Code. **~641★**. | `[verified]` (2026-10-01) |
| `aipoch/medical-research-skills` | Referenced as the validation base for a downstream project (`WHYYS88/skill-compat-layer`), so it appears to exist. | `[unverified]` (indirect evidence only, 2026-10-01) |
| `Aperivue / medsci-skills` | Cited in an earlier draft of our ideas inventory. No matching repository found on GitHub search. | `[unverified]` — could not confirm (2026-10-01) |

## Related ecosystem (context, not CF-specific)

These appeared in the same search space and are useful as design references for how medical/scientific agent skills are structured:

- [`LeoYeAI/openclaw-master-skills`](https://github.com/LeoYeAI/openclaw-master-skills) — curated collection of 1,200+ OpenClaw skills (~2.1k★). A good example of the "catalog of skills" pattern.
- [`runkids/skillshare`](https://github.com/runkids/skillshare) — syncs skills across CLI agent tools (Codex, Claude Code, OpenClaw). Relevant if we want cross-tool distribution of our skills.
- [`bio-xyz/BioAgents`](https://github.com/bio-xyz/BioAgents) (archived) — multi-agent "AI scientist" for biological deep-research; literature analysis + reasoning. Pattern reference for multi-agent research loops.
- [`alvinreal/awesome-openclaw`](https://github.com/alvinreal/awesome-openclaw) — curated list of OpenClaw resources, plugins, and memory systems.

## Specific claims flagged as unverified

These appeared in an earlier draft of our ideas inventory and **could not be confirmed**. They are listed so they can be verified or removed, not as facts:

- `[unverified]` "166 skills" / "869 skills" — specific skill counts; no source found.
- `[unverified]` "250,000 scientists" — scale claim; no source found.
- `[unverified]` `cisco-ai-skill-scanner` — a named tool; no repo or doc found under that name.

> To resolve any of these: find the primary source (repo, release page, or paper), then move it to the table above with an access date and change the label to `[verified]`.

## Cystic fibrosis + AI: what to survey next

This section is intentionally a **work list**, not a set of claims. Concrete, citable CF+AI projects should be added here as they're found and verified.

- [ ] Academic CF+ML work (FEV1 prediction, exome/genomics, sputum imaging) — search PubMed, add with citation + access date.
- [x] CFF's strategic direction — summarized from official CFF pages in [`../cff-goals.md`](../cff-goals.md); year-specific 2026/2027 targets remain unverified from the reviewed public pages.
- [x] Public CF/CFTR code repositories — curated in [`cf-repositories.md`](cf-repositories.md); this is a discovery list, not an endorsement or reuse approval.
- [ ] CF Foundation's AI / data initiatives — continue with official pages only; distinguish published work from inferred opportunity.
- [ ] Open datasets: CFTR genomics, CF registry aggregates (public, de-identified only).
- [ ] Existing MCP servers for biomedical literature / gene databases (e.g., NCBI, Ensembl, UniProt) — note license and rate limits.

**Safety guardrail:** any reference added here must be a **public, de-identified** source. Never link to or store patient-identifiable data.
