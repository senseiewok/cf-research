# Compliance review: public repos, third-party reports, and agents

Status: review, written 2026-10-05. Scope: what the lab does with other people's documents and data, what its public repos may contain, and how agents are used on them. It is a self-review by the agent that did the work, from the catalog and from what the documents themselves say. It is not legal advice; the points marked **Human** need a person's decision, and a lawyer's if the lab ever wants certainty.

## What the lab holds and does

| Thing | What happens | Where it lives |
| --- | --- | --- |
| Registry reports (CFF, UK, Australia, Ireland, Canada, ECFSPR) | The maintainer downloaded copies by hand. Agents read them, extract facts and write notes in the lab's own words | `sources/downloads/` (git-ignored, never committed) |
| CFTR-France variant lists | Public downloads; aggregate counts per version were computed locally | the same folder; only the counts and our words enter the repo |
| Page text and quotations used for extraction | Plain text made from the reports, checked quote by quote | outside the repos, in the untracked Files folder |

## What each source's own text says about reuse

| Source | What it says (read on 2026-10-05) | What the lab does |
| --- | --- | --- |
| CF Foundation (CFFPR) | A "Figure Permissions" paragraph tells readers to ask the Registry team before using its charts or data | Cite and link; no charts; paraphrase |
| UK Cystic Fibrosis Registry | Content may not be used or reproduced in publications without the Trust's permission | Cite and link; no tables; a 10-word limit on shared text |
| Australia (ACFDR) | No licence; enquiries about use or reproduction go to the registry | Cite and link; ask before reproducing tables |
| Ireland (CFRI) | The landing page lists CC BY 4.0; the PDF itself carries no licence and the Zenodo record was not read | Cite; reuse with attribution only after the licence is confirmed |
| Canada (CCFR) | No statement on reuse | Terms unread: cite and link only |
| ECFSPR | No terms recorded in the catalog | Terms unread: cite and link only |
| CFTR-France | Data for educational purposes only, read in the patient's clinical context, no warranty; no licence found | Cite and link; no re-hosting; never present a class as a diagnosis |

## Findings

| # | Finding | What was done | Still open |
| --- | --- | --- | --- |
| 1 | Notes in public repos may only paraphrase and quote briefly (`security-browsing` section 5), but nothing measured it | `tools/sources/check_source_overlap.py` measures the longest run of words a draft shares with any source (default 25 words, 10 for the UK report), with a test that failed under seven deliberate breakages | Run it before committing any note built from a report |
| 2 | The published ECFSPR skill had a 26-word run copied from the 2020 report | Shortened to a paraphrase | |
| 3 | Three catalog notes quoted more than the limits allow (two from the UK report, one from the Australian report) | Paraphrased; the catalog and its generated copy were resynced | |
| 4 | The ECFSPR skill says CC0 in its header but contains short third-party quotations | Added a licence note: the dedication covers the lab's wording only | **Human**: same wording for the other skills if you agree |
| 5 | The reports were read by cloud models (Fable, through Claude Code on the maintainer's account) | Recorded here. The pages contain no personal data, secrets or patient records. What the reviewers read was the extracted text of the reports and the CFTR-France summaries. None of the terms read forbids this, but the UK, ECFSPR and Canadian website terms were not read | **Human**: confirm the account's data-retention and training settings; decide whether sources with unread terms go only to local models |
| 6 | The reviewers were told to use only read-only local tools. The harness did not enforce it; the transcripts show 335 file reads, searches and globs, two directory listings in a shell, and no web tool | Noted as a process gap (ledger E66) | Use an agent type with a tool list when the restriction matters |
| 7 | Registry reports are aggregate; CFTR-France lists variants, not people. Small counts are suppressed by the registries themselves | Nothing to remove. Rule written down: no attempt to rebuild individuals, and the CFTR-France individual database (behind a login) is never touched | |
| 8 | A class in these lists or reports is not a diagnosis, and a registry share is not a person's chance | The new skills say so, and cf-research-context carries the short rule | |
| 9 | robots.txt groups that name AI agents (for example one that disallows `Claude`) did not apply to the evidence client | Decided by the maintainer on 2026-10-05 and done: a fixed list of AI-agent tokens now applies to every `fetch` source (30 new tests, 7 mutants caught); a group naming our own token still decides alone | The list is short on purpose; extend it with a test after checking a token in the operator's own documentation |
| 10a | `fetch_sources.py` downloaded with the library's default user agent and no run-time robots check, relying on the catalog's `robots: allow` field; the registry PDFs it fetched earlier this session went through that code | It now uses the evidence skill's robots rule and conduct (honest user agent, https only, same-host redirects, no retry, a pause, a size cap), with 23 tests and 12 mutants caught; `--audit-robots` re-checked all 10 `fetch` entries on 5 hosts under the stricter rule: all allowed | The earlier requests themselves cannot be re-run under the new conduct |
| 10b | Public notes quoted source sentences, mostly terms-of-use statements in the catalog | Decided by the maintainer: do not quote. 24 catalog entries, the ECFSPR skill and three older notes were paraphrased; `check_source_overlap.py` now also flags any quotation over 8 words | The older notes keep quotations of the lab's own wording |
| 10 | The licences: cf-lab is MIT, cf-research CC BY 4.0, cf-skills CC0, with T-0001 still open for the control repo | Third-party quotations are not covered by any of them, and the licence note in finding 4 says so | **Human**: T-0001; copyright in text an agent wrote may not exist in some jurisdictions, so the lab grants only what it holds |
| 11 | Disclosure of agent work | Commits carry an agent co-author line, and provenance notes name the model that read each source and what a script verified | |

## Not engaged on this evidence

HIPAA and GDPR need personal data; the sources are aggregate, and the individual-level data of the CFTR-France database was never opened. The lab's tools are not a medical device because they diagnose nothing and the notes say so. The EU AI Act's duties for AI systems do not obviously reach a documentation repo, and nothing here was checked against them in depth.

## Rules this review adds

1. Before a note built from a report is committed, run `check_source_overlap.py` on it.
2. Do not quote a source's sentences in a public note: paraphrase and cite report, year and page. Titles, names, labels and field names are short enough to quote. Never a table, a chart's numbers wholesale, or a run of figures. `check_source_overlap.py` enforces this.
3. Terms unread means cite and link only, and no more text than needed goes to a cloud model.
4. A restriction given to a subagent in words is checked against its transcript before the output is used.
5. Fetching follows the evidence skill's network rules, robots.txt groups that name AI agents included. Safe web practice is the first rule, and skills that fetch fork from the evidence skill's `NETWORK-RULES.md`, never from a copy.
6. Every number in a note built from reports is in its evidence or listed with its arithmetic (`check_numbers.py`).
