# Anti-hallucination strategy: what the lab's own record says, and what to build

Status: proposed, written 2026-10-08. The maintainer decides what to adopt. Nothing here changes a rule or a skill; it proposes them. Two analyses fed it: a hand classification of the lab's error ledger (138 rows, one rater, who is also an AI agent, so the counts are lower bounds) and a short literature check (abstracts only, read 2026-10-08). Where a number was re-read by a second agent it says so.

## 1. What the lab's own record says

The ledger (`ledger/`, 138 rows over two sessions) records errors that agents made or that a check caught.

| Kind of error | Rows | First caught by |
| --- | --- | --- |
| Tooling, escaping, environment | 39 | A script in 27 of 39 |
| A check that was vacuous or failed open | 24 | A script 10, the agent 11, another model 3 |
| A wrong check or spec | 19 | Mostly scripts and the agent |
| **A claim stronger than its source** | 9 | **No script.** Agent 3, another model 6 |
| **An invented number, fact or date** | 7 | **No script.** The agent 7 |
| Stale fact (5), confidently wrong model review (4), wrong attachment (3), absence from a partial view (3), premature "fixed" (2) | 17 | A script 1; the agent 10; another model 4; a sub-agent 2 |

Counting the seven claim-type rows together (33 rows): **a script caught 1, the agent caught 20 itself, another model caught 10, a sub-agent 2, and a person caught none first.** The 20 self-catches are discipline, not a mechanism.

Five patterns predict a false claim: numbers or dates typed from memory; working from a summary instead of the primary text; widening words such as "only", "same", "agrees", "stronger"; a verbatim quote treated as proof of the sentence beside it; and a model checking its own claim (two agreeing samples from the same kind of model did not help). Entry ids for each pattern are in the analysis. Five rows were caught only after commit, merge or publication.

## 2. What published work adds (abstracts read 2026-10-08)

| Finding | Source | Read by |
| --- | --- | --- |
| Models "struggle to self-correct" without external feedback and sometimes get worse | arXiv 2310.01798 | one agent, abstract |
| Verifying with questions answered independently of the draft reduces hallucination (Chain-of-Verification) | arXiv 2309.11495 | one agent, abstract |
| Splitting text into atomic facts and checking each against a source is a workable method (FActScore; SAFE agreed with crowd annotators 72% of the time on about 16,000 facts) | arXiv 2305.14251, 2403.18802 | one agent, abstract |
| In generative search engines only 51.5% of sentences were fully supported by their citations and 74.5% of citations supported their sentence | arXiv 2304.09848 | two separate reads, abstract |
| ChatGPT fabricated 55% (GPT-3.5) and 18% (GPT-4) of citations; of the real ones, 43% and 24% had substantive errors | PMC10484980 (2023) | two separate reads, abstract |
| Models guess rather than say "not sure" because training and benchmarks reward guessing | arXiv 2509.04664 | one agent, abstract |

Limits: abstracts only, and results on other tasks and older models; they point to where to look, not to the lab's own error rates.

## 3. Principles

1. **Ground in primary text, and say which kind you used.** A search summary or a news item is a lead, not evidence. Every answer carries a tier (section 4).
2. **A claim is a sentence plus an exact quote, a source id and a date.** No quote, no publication. A quote supports only what it says: its years, no cause unless the source gives one.
3. **Widening words need a scope record.** "Only", "all", "never", "same", "first", "no longer", "new", "agrees", "stronger" and absence claims must name what was searched and where. This is the class no script catches today.
4. **Numbers come from a source line or a command.** Arithmetic is written out. (In place; extend it to every note and page.)
5. **Independent verification.** The verifier sees the question and the source, never the draft or the author's verdict; it is a different model family when a model is used; a script runs first; every flag is checked on the source; a majority is never a verdict.
6. **"Not stated" is a correct answer.** Score it above a wrong answer in every benchmark and prompt.
7. **Checks fail closed.** A check that skips what it cannot apply, loops over nothing, or reads empty results as an answer must fail. Every new checker ships with a negative control.
8. **Control drift.** One source of truth per fact, generated views, counts from commands, dated notes. A sub-agent's report is a lead, re-checked on the source; a rule given to an agent only in words is checked against its transcript.

## 4. Answer tiers (for chat, notes and pages)

| Tier | Meaning | May be used for |
| --- | --- | --- |
| T0 | From a summary, search snippet or memory. Unverified | Orientation only; label it |
| T1 | Quoted from a primary source, with id and date | Notes and pull requests |
| T2 | T1 plus a script check and a different model's check, flags confirmed on the source | Public claims |
| T3 | T2 plus a person's review | Anything that could affect a person's care, any scientific copy on the site |

A patient-related question needs T1 at least before an answer, and says so. This session's CF answers began at T0: a search summary said gut changes were "linked to growth" and the primary paper said the opposite.

## 5. What exists, and the gap

| Failure | Defence in place | Gap |
| --- | --- | --- |
| Number not in the evidence | `check_numbers.py` | Proves a number exists, not that it sits on the right claim |
| Copied text | `check_source_overlap.py` | Copying, not truth |
| Citation that does not resolve | Evidence loop (single lookups) | No scan of a whole note |
| Claim stronger than its quote | A human or another model, sometimes | **No mechanism** |
| Absence claim from a partial view | Agent discipline | **No mechanism** |
| Vacuous or fail-open check | Review | No rule that skipped means failed |
| Model review confidently wrong | Reading the source | Keep: flags are leads |

## 6. Build order (ranked by what the ledger says matters)

1. **Claim ledger checker** (`tools/claims/`, cf-research): reads a claims file (id, sentence, source id, exact quote, scope note) and the evidence text; fails when a quote is not an exact substring, a number in the sentence is not in the quote or a listed computation, a widening word has no scope note, a claim has no quote. It cannot prove entailment (the quote still has to support the sentence), but it removes the cheap classes and forces the scope record. Negative controls from ledger cases (E59 a number from memory; E73 a sentence past its quote; an absence claim with no scope).
2. **Cite-check** (`cf-evidence-loop`): scan a note for DOIs, PMIDs and trial ids; resolve each through the existing providers; flag a missing one or a title that does not match. Same hosts, no new dependency. Board T-0010's frozen claim set is the matching regression suite.
3. **Frozen claim set with known-false variants** (board T-0010): rerun after any model or prompt change; scores "not stated" above wrong.
4. **Fail-closed audit of the existing checks**: the mutation harness skips mutants it cannot apply and still passes (ledger E89, E111, E137); empty search results read as answers (E8, E121); checks that loop over nothing (E133). A skipped item becomes a failure unless allowed by name.
5. **Verifier-blind packets**: the claim-check packet carries the claim and the source only. Already the lab's rule for some packets; make it the template.
6. **New grounding sources**: RxNorm/RxNav (terms read: free, 20 requests per second per IP, attribution required) for drug names; PubTator3 and the EBI ontology lookup for gene, variant, disease and phenotype normalisation (terms unread); OpenAlex and MyVariant are catalogued but have no command. Each needs a catalog entry first and a person to admit it.
7. **Canary strings and red-team packets**: no source read; treat as a hypothesis.

Not solved by any of this: whether a true quote supports a sentence in the reader's sense (a person or a different model, then the page), anything nobody tests, and errors that were never caught (the ledger cannot count them).

## 7. Decisions for the maintainer

- Adopt the tiers (section 4) as the lab's wording for what an answer rests on?
- Adopt principle 3 (a scope record for widening words) as a rule in `cf-research-context`? It is human-maintained policy, so it is a proposal here, not an edit.
- Build in the order of section 6? Items 1 and 2 are being built now as separate pull requests (they add tools only and change no rule); the rest wait for a decision.
