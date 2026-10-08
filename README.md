# Sensei Ewok Research

Open-source research notes from Sensei Ewok's CF Lab. This is a
knowledge base and an ideas inventory, kept openly so a stranger can check
it. Alongside the notes it holds small, single-purpose tools under `tools/`, built
so the CF community can run and check them. It is not a clinical product.

```text
⠀⠀⠀⠀⠀⡀⠄⠂⠂⠁⠁⠛⠛⠶⢦⣄⡀⠀⠀⠀⠀⠀⠀⠀⡀⠄⠂⠂⠁⠁⠛⠛⠶⢦⣄⡀⠀⠀⠀⠀⠀⠀⠀⡀⠄⠂⠂⠁⠁⠛⠛⠶⢦⣄⡀
⠀⠀⠀⠄⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⠶⣄⠀⠀⠀⠄⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⠶⣄⠀⠀⠀⠄⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⠶⣄
⠀⢦⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⢦⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⢦⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⢦
⠀⠀⠈⠳⢦⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠂⠁⠀⠀⠈⠳⢦⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠂⠁⠀⠀⠈⠳⢦⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠂⠁
⠀⠀⠀⠀⠀⠉⠛⠶⢦⣤⣤⡄⠄⠂⠂⠁⠀⠀⠀⠀⠀⠀⠀⠀⠉⠛⠆⠀⠀⠄⠠⠀⠒⢆⠁⠀⠀⠀⠀⠀⠀⠀⠀⠉⠛⠶⢦⣤⣤⡄⠄⠂⠂⠁
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠑⠢⠄⣀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠂⠒⠤⠄⠤⣀⡀⣀⡀
```

*A double helix drawn in Braille dots twists gently across an empty field; one strand has a small gap near the middle, and a thin thread curves up toward it from the lower right.*

> *Into the gap we work a thread of questions, slowly, and check each one before the next.*

**Research notes, a catalog of 53 sources and what each permits, a 3D model of the CFTR protein from four published structures, and a ledger of our own mistakes.** Research, not medical advice.

## What this is

This repository is the lab's notebook, kept where anyone can read it. The notes are plain Markdown and need no toolchain. Every claim in them that a stranger could check carries a label: `[verified]` with its source and the date we checked, `[unverified]`, or `[hypothesis]`, so you can see at a glance what we know, what we think and what we have not confirmed. [How to read claims](#how-to-read-claims) explains the three.

Beside the notes are the things the notes rest on. The **source catalog**, [`sources/catalog.yaml`](sources/catalog.yaml), records every external source we rely on and what its publisher permits. On 2026-10-05 it held 53 sources, counted from the file: 15 a tool may fetch after reading the site's `robots.txt`, 14 APIs with a recorded usage policy, 20 a person must download by hand, 3 refused by their publisher and 1 that needs a written request. Our tools may use nothing the catalog does not allow, and no source document is committed here; they are linked and checksummed, never redistributed. Three **reading guides** under `.claude/skills/` say how to read the CFF and ECFS registry reports, and genotype and variant figures in registry reports and variant lists, before quoting a number from them. No registry scientist has reviewed the guides yet, so they are reading aids, and we compare no figure from one registry with another.

The **CFTR model**, [`tools/cftr-model/`](tools/cftr-model/README.md), is an interactive 3D drawing of the human CFTR protein built from four published structures (6MSM, 5UAK, 8EIQ and 8EJ1, from the Protein Data Bank), numbered as in UniProt P13569. It draws one point per amino acid each structure places and says plainly what it leaves out: no pore, no gate, no medicine, no motion, and nothing about what any person with CF can expect. Residue 508 is marked where a structure has it and shown as a gap where the file says it is absent, as it is in two of the four. The model is a draft, in review; the website's [CFTR model](https://senseiewok.ai/cf/cftr/) page is the same model with the same tables. The picture at the top of this page echoes that gap, but it is an emblem of the work of understanding, not a diagram of a gene or a protein, and it shows no treatment.

The **ledger**, [`ledger/`](ledger/), logs the real mistakes caught in our working sessions, what caught each one and a cheap guard against the next: 138 entries in two session logs on 2026-10-08, by `ledger/tally.py`. A single entry is an observation, not a lesson. [`proposals/`](proposals/) holds changes to sibling repositories that an agent could not or should not apply directly, each with its evidence and a verification step.

What this is not: a medical device, clinical guidance, or a place for patient data. We do not interpret anyone's genotype, recommend or dose any medicine, or predict any person's outcome. Clinical decisions belong with licensed care teams.

## Try it

The notes need nothing: open any `.md` file. The tools need Python; each has its own README and dependency list under `tools/`.

```sh
python tools/cftr-model/serve.py                 # the CFTR model at http://127.0.0.1:8080/ , this computer only
python tools/sources/fetch_sources.py --manual   # the sources a person must download by hand; prints a list, fetches nothing
cd tools/sources && python -m unittest -v        # the source tools' tests, no network; needs PyYAML (see tools/sources/requirements.txt)
```

The model needs a browser with WebGL2; without it the page still reads, because every fact is also written out as text and tables. Rebuilding the model's data from the four structure files, and the browser tests that need Playwright, are in [`tools/cftr-model/README.md`](tools/cftr-model/README.md); the fetcher's options and what it will not do are in [`tools/sources/README.md`](tools/sources/README.md).

## What's inside

The [Layout](#layout) table below lists every path with what it is and its status. In short: the notes at the top level and under `landscape/`; the catalog and its rules under `sources/`; the tools under `tools/` (`cftr-model`, `sources`, `skill_trial`, and a pointer for the evidence loop, which moved to the skills repository); the three reading guides under `.claude/skills/`; `proposals/`; `ledger/`.

## How we check it

- **Labels on claims.** `[verified]` carries a source and a date; anything else is not to be cited as fact.
- **Scripts that fail on invented text.** `tools/sources/check_numbers.py` flags every number in a draft that its evidence does not contain, and arithmetic must be listed with its sum (`564  # 127 + 437`) so a person can check it. `check_source_overlap.py` measures the longest run of words a draft shares with the report it was written from, fails above a limit (stricter for a registry that asks permission before its content is reproduced) and flags quotations longer than eight words. A pass means the numbers exist in the evidence, not that each sits on the right claim; a second model still reads each claim beside its evidence.
- **The model's numbers come from files.** `tools/cftr-model/build_traces.py` reads the structure files and UniProt, compares every placed residue with UniProt's sequence and refuses to continue if one disagrees without its file saying so. Its tests check that the fit recovers a known rotation, that the generated page blocks match the data and, in a real browser, that every answer a research tool gives equals a value recomputed from the data. Tests that can never fail prove nothing, so the guards were also checked by breaking the code on purpose; the record is in the tool's README and in [`proposals/2026-10-05-cftr-model-review.md`](proposals/2026-10-05-cftr-model-review.md).
- **Reproducible sources.** `sources/manifest.json` records the SHA-256 of the exact bytes an extraction ran against, so a stranger can download the same report from its publisher, confirm the hash and reproduce the result without us hosting a copy.
- **A ledger, not a scoreboard.** Every caught mistake is recorded with what caught it. The website's [Evidence](https://senseiewok.ai/evidence/) page shows the catalog counts, one evidence record and three ledger entries exactly as recorded.

## Working with the other lab repositories

This repository is one of three that sit side by side; the control repository, [`cf-lab`](https://github.com/senseiewok/cf-lab), holds the shared agent rules. Start agent sessions from the `cf-lab` folder where you can, and add this one as a working directory: its section [Working across the sibling repos](https://github.com/senseiewok/cf-lab/blob/main/AGENTS.md#working-across-the-sibling-repos) says how, and what does and does not load.

If you do start a session here, `.claude/settings.json` keeps the same secret-file deny rules as `cf-lab`, so the agent may not read `.env` files, keys or credential folders. A change that spans repositories is one pull request in each.

## To our CF community

To people living with cystic fibrosis, families, caregivers and researchers: this is the notebook of a small, independent lab, and you are welcome to read over our shoulder. What is here is for examining evidence: notes that say how sure we are, a record of what each source permits, a model of a protein that shows only what four published files show, and a list of the mistakes we have made, so you can judge for yourself how carefully we work. None of it is medicine or medical advice; we do not interpret anyone's genotype or say what any person can expect, and clinical decisions belong with a care team. The sixty-five roses on the lab's [CF story](https://senseiewok.ai/cf/) page are our tribute, independent of the Cystic Fibrosis Foundation and not an emblem. You deserve care, dignity and room for ordinary life, and you owe no one an inspiring story. One concrete ask: the CFTR model is a draft, and we would like a person with CF, or someone who cares for one, to read its wording before it leaves draft. If a sentence there, or anywhere here, is unclear or lands wrong, open an issue in this repository and quote it; we will change it and say what changed.

With care,

Sensei Ewok

```text
                        (   )
                     (  * * *  )
                   (  (       )  )
                (  (   (  *  )   )  )
               (  (  (  (    )  )  )  )
              (  (  (  (  *  )  )  )  )
             (  (  (  ( ( ( ) )  )  )  )
             (  (  ( ( ( ( ) )  )  )  )
              (  (  ( (   )  )  )  )
               (  (  (  )  )  )
                (  (   )  )
                 (  (  )
                   ( )
                    (
                    (
                    (
```

> *You will cross the desert to find what you already carry. That is not a
> disappointment. That is the whole point of crossing.*

The notes below are that road: what we know, what we think, and what we have not
verified yet, in that order. Read them the way you'd read a map, not a promise.

> **Not medical advice.** Nothing here diagnoses, treats, or recommends dosing. No patient data (PHI) belongs in this repo. Anything touching a real patient goes through care teams. See the [disclaimer](#disclaimer).

## Layout

| Path | What it is | Status |
| --- | --- | --- |
| [`cystic-fibrosis.md`](cystic-fibrosis.md) | Disease knowledge base — pathophysiology, CFTR, clinical landscape | v1 |
| [`cff-goals.md`](cff-goals.md) | CFF's 2026–2030 plan, community input, and what is or is not public for 2026–2027 | research note, checked 2026-10-03 |
| [`cf-ai-ideas.md`](cf-ai-ideas.md) | AI agent / skill / MCP ideas for CF research, organized by tier and ROI | draft |
| [`landscape/medical-ai-landscape.md`](landscape/medical-ai-landscape.md) | Verified external projects, repos, and references in medical AI / agent skills | draft |
| [`landscape/cf-repositories.md`](landscape/cf-repositories.md) | Curated GitHub CF/CFTR research and software repositories, with license and validation cautions | snapshot, checked 2026-10-03 |
| [`landscape/registry-data-access.md`](landscape/registry-data-access.md) | What the CF registry publishers do and do not permit programmatically, where to download their reports, and a draft permissions email | tested 2026-10-04 |
| [`sources/catalog.yaml`](sources/catalog.yaml) | Every external source we rely on, with its access permission (`fetch`, `api`, `manual`, `request`, `forbidden`), URL, and verification status. Agents read this instead of guessing | v1, 2026-10-04 |
| [`sources/README.md`](sources/README.md) | The rules an agent must follow when reading the catalog, and why no source PDFs are committed | v1 |
| [`landscape/antibiotic-stewardship-cf.md`](landscape/antibiotic-stewardship-cf.md) | A short argument for antibiotic stewardship in CF, its five claims each checked against the primary papers, with exact quotes, links and limits | research note, read 2026-10-07 |
| [`landscape/cftr-structures.md`](landscape/cftr-structures.md) | The four published structures of human CFTR behind the model, how they were fitted onto one another, what each file declares, and what the model does not show | research note, checked 2026-10-05 |
| [`tools/cftr-model/`](tools/cftr-model/README.md) | An interactive 3D model of CFTR from those structures (WebGL2, no libraries), research tools (look up a residue, measure a distance, layers from UniProt, export), the build that ties every number to its file, and its tests. Draft, in review | draft, in review |
| [`tools/sources/`](tools/sources/README.md) | Fetches what the catalog permits, prints browser instructions for the rest, writes a SHA-256 manifest for reproducibility | working, smoke-tested 2026-10-04 |
| [`tools/skill_trial/`](tools/skill_trial/) | With-and-without trials of agent guidance on the local worker, and the checkers that judge them (page basics, ASCII art, skill discovery). Used by [`proposals/2026-10-07-public-skills-trial.md`](proposals/2026-10-07-public-skills-trial.md) | working, 2026-10-07 |
| [`.claude/skills/ecfspr-report-reading/`](.claude/skills/ecfspr-report-reading/SKILL.md) | Agent skill: how to read an ECFSPR annual report (2020 to 2024) before extracting numbers. Every quotation and layout claim was checked against the report text; see its provenance note | draft, content review pending |
| [`.claude/skills/cffpr-report-reading/`](.claude/skills/cffpr-report-reading/SKILL.md) | Agent skill: how to read a CFFPR Annual Data Report (2020 to 2024) before extracting numbers: who is counted, how a value is made, what changed between years and what the reports warn about | draft, content review pending |
| [`.claude/skills/cf-genetics-reading/`](.claude/skills/cf-genetics-reading/SKILL.md) | Agent skill: how to read genotype and CFTR-variant figures in six registries' reports and the CFTR-France variant lists before extracting or comparing a number | draft, content review pending |
| [`tools/evidence/`](tools/evidence/README.md) | Pointer: the evidence loop moved to the portable skill `cf-skills/.claude/skills/cf-evidence-loop`. Run it here with `EVIDENCE_CATALOG` set to `sources/catalog.yaml` | moved 2026-10-04 |
| [`.claude/settings.json`](.claude/settings.json) | Claude Code deny rules for a session started here: no reading `.env` files, keys or credentials. The same list as `cf-lab` | v1, 2026-10-06 |
| [`proposals/`](proposals/) | Proposed changes to sibling repos that an agent could not or should not apply directly, with evidence and a verification step for each | — |
| [`ledger/`](ledger/) | Logs of real errors from working sessions: what happened, what caught it, and a cheap guard. A single entry is an observation, not a lesson | two session logs, 2026-10-04 and 2026-10-07 |

## How to read claims

Every factual claim that could be checked by a stranger is labelled with one of:

- **`[verified]`** — confirmed against a live source (repo, citation, or official page) on the date written. The source + date are given inline.
- **`[unverified]`** — plausible but not independently confirmed. Treated as a hypothesis. Do not cite these as fact.
- **`[hypothesis]`** — our own idea or inference, not a claim about the external world.

If you add or edit a claim, update its label and, for `[verified]`, keep the source link and the date you checked it.

## Contributing conventions

- Markdown, sentence-case headings, short paragraphs, tables over prose where they aid scanning.
- Prefer editing an existing file over creating a new one; add a new file only when there's no existing home.
- Notes stay readable with no toolchain. A stranger should be able to read every `.md` here without installing anything.
- Code lives under `tools/`, one directory per tool, each with its own dependency list and README. Nothing outside `tools/` needs dependencies.
- Never commit secrets, patient data, or third-party report PDFs. Source documents are linked and checksummed, not redistributed — see [`landscape/registry-data-access.md`](landscape/registry-data-access.md).
- Reference external work by stable URL + access date, not by star count alone (star counts change).

## Disclaimer

This repository contains research notes, ideas, and references. It is not a medical device, not clinical guidance, and does not contain or request patient-identifiable data. Clinical decisions belong with licensed care teams.

## Licenses

Notes and documentation are licensed under [CC BY 4.0](LICENSE). Code under `tools/` is MIT ([tools/LICENSE](tools/LICENSE)). The three agent skills under `.claude/skills/` each declare CC0-1.0 in their frontmatter. Source documents published by registries are not covered by these licenses and are not redistributed here; short quotations are attributed to their report and year.
