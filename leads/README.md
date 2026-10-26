# Research leads

A **research lead** is an idea, a question or a source worth following that is not yet a task. This folder is the lab's
public register of them. It sits upstream of the lab's task board: nothing here is planned work, and nothing here is a
finding.

| File | Role |
| --- | --- |
| `leads.json` | The single record. Edit this |
| `LEADS.md` | The readable view, rendered by `tools/leads/render_leads.py`, grouped by status then id. Never edit it by hand |
| `../tools/leads/check_leads.py` | Checks `leads.json` against the rules below and against `sources/catalog.yaml` |

```sh
python tools/leads/render_leads.py            # after editing leads.json
python tools/leads/check_leads.py             # exit 0 pass, 1 a rule failed, 2 usage or input problem
python tools/leads/render_leads.py --check    # exit 1 when LEADS.md is stale
python -m unittest discover -s tools/leads -v # the tools' tests
```

Standard library only; no network.

## From a lead to work

1. **Lead.** Anyone, a person or a model, adds a lead here with `status: idea`.
2. **Scoping note.** When a lead looks worth doing, someone writes a short note in `proposals/` (named
   `YYYY-MM-DD-short-name.md`) with the question, the data and what the catalog permits, the guard and a check. The lead
   moves to `scoping` and links the note.
3. **Board row.** Someone writes a `proposed` row on the lab's task board (`tasks/board.json` in the control repository)
   and links it here. A model may write a `proposed` row; **only a person sets a row `ready`**.
4. The lead follows the work: `in-progress` and `done` each need a board row or a pull request in `links`.

`data-available` means the data the lead needs is already in hand (for the trial atlas, in a snapshot we have), not that
anyone has started.

## Fields

| Field | What it holds |
| --- | --- |
| `id` | `L-NNNN`, never reused |
| `title` | A short name |
| `question` | One sentence |
| `why_it_matters` | One sentence, for a CF researcher or clinician |
| `status` | `idea`, `scoping`, `data-available`, `in-progress`, `done`, `dropped` |
| `sources` | A list of `{type, ref, tier, read_on, note}`. `type` is `registry`, `literature`, `data`, `council`, `conference-private` or `run`. `tier` is T0 to T3; T1 or higher needs `read_on`. `note` is short: never a long quote (at most 25 words in quotation marks). Optional `read_via` says how the source was read (see below) |
| `data_needed` | `{catalog_ids, access, already_have, fields_missing}`, and optionally `not_in_catalog` (plain names of sources the catalog does not list). Every id must be in `sources/catalog.yaml`, and `access` must be the catalog's value (`api`, `fetch`, `manual`, `request`, `forbidden`); when the ids differ, `access` is a mapping of id to value. With no catalog ids, `access` is `not in catalog` |
| `effort` | `{size: S, M or L, estimate: true}`. Always an estimate |
| `guard` | `{ladder, note}`: how every public number or sentence from the lead would be grounded. `ladder` uses `computed`, `registry-statement-verbatim`, `open-access-quote`, `lab-position` |
| `must_not_claim` | The over-claims this lead invites. Never empty |
| `needs_person` | `clinician-reader`, `registry-scientist`, `maintainer-approval`, or `none` alone |
| `links` | `{board: [T-NNNN], proposals: [proposals/...md], prs: [repo#N]}` |
| `created`, `updated` | `YYYY-MM-DD` |
| `created_by` | `{kind: human or model, model_id}`; a model lead carries its `model_id` |
| `visibility` | Always `public` in this file |
| `notes` | Optional, short |

## Tiers: say how sure we are

Tiers follow the lab's grounding rules: **T0** a summary, snippet or memory, unverified; **T1** a quote from a primary
source with its id and date; **T2** T1 plus a script check and a different model's check; **T3** T2 plus a person's
review. A lead states each source's tier honestly. **A T0 lead stays labelled T0 until someone reads a primary source**,
and then the source gains `read_on` and its new tier; never upgrade a tier on the strength of a summary. A number a
model read in an abstract is T0 here, and the `note` says so. A verified quote proves only that the words exist; a
person or a different model still reads it against the sentence it supports. `LEADS.md` shows each lead's best source
tier; every other source keeps its own.

### How a source was read

`read_via` records it, and the checker holds the tier to it:

| `read_via` | Meaning | Tier rule |
| --- | --- | --- |
| `direct` | The page read as it is | Any |
| `local-text-extract` | A PDF or file turned into text on this machine and read | Any |
| `lab-run` | Output of the lab's own script over a dated copy of the data; the note says whether it was recomputed independently | Any |
| `page-summarising-tool` | A page read through a fetch tool that summarises pages, so the quote is as the tool returned it | At most T1, and the note must contain: `read through a page-summarising tool; quote not checked against the raw page` |
| `search-snippet`, `memory` | A search result snippet or something remembered | T0 |

A source whose note carries that sentence must have `read_via: page-summarising-tool`. Before a tool-read quote goes
into public text, someone reads it on the raw page.

A source under a non-commercial licence (CC BY-NC and its variants, named in `ref` or `note`) is paraphrased and
cited, never quoted, in this file: such a quote is private evidence only. The checker refuses quotation marks in the
note of such a source.

### Sources no tool may use automatically

Some sources a lead might want are closed to automation by the catalog. **CFTR2** is the example that comes up most
for genotype leads: `sources/catalog.yaml` records it as `manual`, link-only and not redistributable (T1 to the
catalog's notes, read on 2026-10-05; its terms were not re-read for this register). So no lead plans an automated
CFTR2 lookup; a lead may link to it for a person to read. The same holds for any source whose catalog `access` is
`manual`, `request` or `forbidden`.

## Public register, private folder

This file is public. A lead that **cannot** be public does not go here at all:

- third-party material: anything from a conference poster or talk, an unpublished manuscript, a private message;
- anything that identifies a person or an institution through what they presented or said;
- anything else private: patient data, site details, keys, transcripts, email addresses.

Such leads stay in the maintainer's private lab-files folder (the untracked `.lab-files` folder inside
`cf-lab`, in its `memory/` area, under the rules in the control repository's `AGENTS.md`, section "Local files and
memory"). The public register may point at one only through a source of type `conference-private` with `ref: "private"`
and an empty note, as L-0019 does: no names, no institutions, no poster contents. The checker refuses any other detail
on such a source, any email-like text, any user folder path, and any trial registration id.

## Adding or updating a lead

1. Take the next free id. Fill every field; write plainly, and say "not stated" rather than guess.
2. For a source the lead needs, find its id in `sources/catalog.yaml` and copy its `access`. If it is not there, list it
   in `not_in_catalog`; using it needs a catalog entry, with its terms read, first.
3. Set `updated` to today when anything changes.
4. Run the renderer, then the checker, and commit both files together.

## Dropping or archiving a lead

Never delete a lead: the record of what was considered is part of the point. Set `status: dropped` and say why in
`notes` (for example: the data is forbidden, a published source already answers it, or it cannot be grounded). A `done`
lead stays as it is, linked to the work that finished it.
