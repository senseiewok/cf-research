# Proposal: a website section for the evidence loop

Status: proposal for the `.ai` repo (board T-0041). Draft copy in `lab-voice` register; every number below resolves to a record the tool produced on 2026-10-04. Not applied; website edits are Elevated tier.

## The name

The site already has "The AI-Loop Field Guide" for how the lab runs models. This is its counterpart for how the lab answers a research question. Options considered:

| Name | For | Against |
| --- | --- | --- |
| **The evidence loop** | Pairs with "AI loop"; says what it is; plain | None serious |
| The researcher's loop | Warm; names the role | Implies a person the lab does not employ |
| The evidence desk | Journalistic, concrete | "Desk" suggests a service you can write to |
| Claim check | Short | Sounds like a fact-checking brand; narrower than the seven layers |

Recommendation: **The evidence loop** for the site section; `evidence` is the package name; the skill is `skills/cf-evidence-loop`.

## Where it goes

After "The AI-Loop Field Guide", before "From the lab notebook". Same card structure as the field guide so the two read as a pair: one explains how the machines work, the other explains what they are pointed at.

## Draft copy

> ### SE-LAB · Vol. 02
> ### The evidence loop
>
> A plain-language guide to how this lab answers a research question, and why every answer comes with its own limits attached.
>
> **1. Seven questions, in order**
>
> A researcher does not start with an answer. They walk a fixed set of questions: What has been published? Has any of it been retracted or corrected? What is being tested in trials, on whom? What is actually approved, where, and for whom? How many people does this concern, and how are they doing? What does the field recommend? And what is the field saying to itself, in grants, preprints, and conference halls, before any of it is published?
>
> We built a small tool that asks those questions of the public record, one source at a time. It does not think. It fetches, copies the exact fields it relied on, stamps the date, and writes down what the answer does *not* establish.
>
> **2. What an answer looks like**
>
> Here is one, unedited, from the day we built it. We asked whether a 2013 paper that underpins the CFTR2 variant database has ever been retracted or corrected.
>
> ```
> status        found
> source        Crossref, with Retraction Watch data
> doi           10.1038/ng.2745
> notices       0   (clean)
> accessed      2026-10-04
> limitations   Crossref holds no update notice for this DOI at access time.
>               Publishers lag and some notices lack the relation, so absence
>               here is strong but not conclusive. Does not assess the quality
>               of the work.
> ```
>
> The last field is the one we care about most. A clean result is good news with a footnote, and the footnote travels with it.
>
> **3. The same tool, pointed at a known case**
>
> Asked about a 1998 paper that was later withdrawn, it returned two notices: a partial correction in 2004 and a full retraction in 2010, each with the notice's own DOI. No judgement, no summary. Two records, two dates, and a pointer to the publisher's words.
>
> **4. What it will not do**
>
> It reaches only sources whose terms we have read and recorded. If a publisher says no to automated access, the tool cannot be made to try; that refusal is a test in the codebase, not a policy in a document. It never presents itself as a browser. It never calls a language model. And it never touches data about an individual person.
>
> **5. Why this is the lab's standard**
>
> Every number on this site, in our notes, and in any video we publish is held to the same shape: source, date, exact field, limitation. When we cannot produce that record, we do not publish the number.
>
> *Research, not medical advice. The evidence loop checks public records; it does not interpret them for any person. Decisions about care belong with you and your care team.*

## Numbers in the copy, and where they resolve

| Claim in copy | Record |
| --- | --- |
| 2013 paper, DOI 10.1038/ng.2745, zero notices | `evidence/ledger-smoke.jsonl` line 1, and `skills/cf-evidence-loop/tests/fixtures/crossref_updates_clean.json` |
| 1998 paper, correction 2004, retraction 2010 | `skills/cf-evidence-loop/tests/fixtures/crossref_updates_retracted.json` |
| "seven questions" | `skills/cf-evidence-loop/references/ARCHITECTURE.md`, layer table |

The withdrawn paper is not named in the copy on purpose. Naming it invites a conversation the lab does not need to host; the point is the tool's behaviour, not the paper.

## Before it ships

- Human read against the `lab-voice` review list, especially "Are the lab's aspirations clearly separate from tested outcomes?" The copy describes a working tool with 17 tests; it should not describe anything the tool does not yet do.
- The `.ai` repo's own `write-like-sensei-ewok` skill was not consulted. If it differs from `lab-voice`, decide which governs.
- Add a Playwright check that the section renders and the disclaimer block is present, alongside the 65-rose count (T-0016).
