# Small AI tools for CF research

**Proposal inventory, refreshed October 6, 2026.** The first version of this page (October 4, 2026) named three starting points and one experiment. Two days of work later, this version says what each became, what we measured along the way, ranks what we could do next, and lists some sparks that need almost no machinery. These are still hypotheses about useful work, not validated community needs or clinical products. Check existing alternatives and seek qualified feedback before choosing a project. Research, not medical advice.

Every count below was read from a file or a command's output on 2026-10-06; the evidence table that lists each one with its source is kept privately. Words such as "observed", "computed", "inferred" and "hypothesis" mark which step each claim is.

## 1. Where we are

### What the October 4 proposals became

| Proposal (October 4) | Status today | Evidence (observed) |
| --- | --- | --- |
| Literature evidence cards | Partly shipped, as a different shape | The `cf-evidence-loop` skill (in the skills repo, version 1.0.0) answers 15 fixed questions against 7 public hosts and returns records, not answers. A record has 7 fields; the code refuses to build one with an empty `limitations`. What did not ship: the study-level card (population, endpoint, result) and the frozen-reference comparison of local work against frontier planning. No paired run against a cloud model has been made |
| Historical claim ledger | Not started as proposed | No sourced timeline exists. What shipped instead: an error ledger of our own mistakes (98 entries, `ledger/tally.py`), and three reading guides under `.claude/skills/`, two of which record how a registry's definitions and layout changed between its 2020 and 2024 reports, while the third records how a variant list's classes changed. That is a ledger of methods, not of claims |
| Research-tool catalog | Partly shipped | `sources/catalog.yaml` holds 53 sources with what each publisher permits: 15 may be fetched after `robots.txt`, 14 are APIs with a recorded policy, 20 must be downloaded by hand, 3 are refused and 1 needs a written request. `landscape/cf-repositories.md` lists 6 repositories with activity and licence metadata and a review note each. A comparison of data-access terms across tools was not built |
| Proposed first experiment | Not run as designed | The experiment asked for frozen fields, independent reference answers and a count of accepted cards. What ran instead were loop measurements on code tasks (below), which say nothing about evidence cards |

### What we learned (measured; small samples; one lab, one machine)

1. **The checks that catch things are mechanical or blind.** Of 98 ledger entries, the first detector was our own re-read 45 times, an interpreter or test 20 times, a script assertion 12 times, a blind review by a different model 7 times and a person twice. Class C (an inference beyond the evidence, or a verifier that could not fail) was the largest class with 43 entries, and the costliest.
2. **A quiet review is not coverage.** The local worker read an 808-line script with 17 confirmed defects four times with a bare packet and found 1, 0, 2 and 2 of them; 9 of the 17 were found by no run. On a 17,000-character note, a fast pass found nothing where a cloud reviewer found 10 points. A finding counts only with a failing input that reproduces (observed).
3. **Cutting a task into verified steps worked where one shot did not.** The same page that failed 0 of 6 attempts in one packet passed 1 of 3 runs when retried whole and 3 of 3 chains as five ordered steps, each with its own verifier (one task; a direction, not a rate).
4. **A data boundary matters more than the model.** Fast mode followed 9 of 16 planted instructions in a plain prompt and 0 of 16 when the text was wrapped as untrusted data and the model had no tools.
5. **Rephrasing a request did not help.** A four-task, three-pass relay pilot gave direct local work 3 of 12 accepted, one preparer plus local 3 of 12 and another 0 of 12, across 60 local calls. Keep the original request; repair from verifier output.

### What is still hypothesis

Who would use any of this, and whether an existing tool already serves them: not known; to be asked of the community. The reading guides have not been read by a registry scientist. The CFTR model page is a draft until a person with CF or a carer has read its wording (board task T-0081). The board holds 97 tasks, 96 of them proposed by a model and 1 by a person, and none has a measured return yet, so every effort estimate on it is still an estimate.

## 2. Ranked by impact

Eight ideas, ordered by expected value for CF researchers, clinicians' educators, families and community writers, per unit of effort. **How we scored (inferred, one author's judgement, not a measurement):** four marks from 1 to 5 each. Reach: how many kinds of reader could use it. Trust gained: how much it lets a stranger check us. Effort: 5 means a weekend, 1 means weeks. Risk: 5 means almost no way to mislead anyone, 1 means a real chance. The total is the sum; ties are broken by reach, then trust. Anyone may re-score it in a pull request; the point is that the scoring is visible.

| Rank | Idea | Reach | Trust | Effort | Risk | Total |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | Claim check before publishing | 5 | 5 | 3 | 3 | 16 |
| 2 | How to read a CF paper, a checklist | 4 | 3 | 5 | 4 | 16 |
| 3 | Replication cards | 3 | 5 | 4 | 4 | 16 |
| 4 | A negative-results shelf | 2 | 4 | 5 | 5 | 16 |
| 5 | An open question board | 3 | 3 | 4 | 5 | 15 |
| 6 | A record behind every number we already published | 2 | 5 | 3 | 5 | 15 |
| 7 | Wider evidence coverage: EU records | 3 | 4 | 2 | 3 | 12 |
| 8 | What changed between registry report years | 3 | 4 | 2 | 2 | 11 |

Each idea below gives the same things: who it could help (a hypothesis), the smallest useful artifact, how we would check it, what is out of scope, what infrastructure and first week it needs, and a stop rule. In every case, a second model's agreement is not a reference answer.

### 1. Claim check before publishing

Someone pastes one sentence they intend to publish and gets back the evidence record that supports it, or an honest "not found".

- **Who (hypothesis):** community writers, educators and the lab itself.
- **Smallest artifact:** a documented workflow that maps a sentence type to one of the 15 existing commands (a DOI goes to `retractions` and `doi`; a brand name to `approval` and `label`; a trial id to `trial`) and prints the record with its `limitations` line first.
- **Check:** a frozen set of sentences with known records, written by hand before any run; count records returned, wrong records, and sentences the tool correctly declines.
- **Out of scope:** judging whether the sentence is true; any sentence about a named person; anything the catalog does not permit.
- **Infrastructure and first week:** a script around the existing tool, no model; write the frozen set and the mapping table.
- **Stop rule:** if more than a handful of frozen sentences produce a confident-looking wrong record, stop and redesign the output so "not found" is unmistakable.

### 2. How to read a CF paper, a checklist

A one-page, plain-language checklist built from our own guardrails: what the denominator is, which year the title page says, whether a number is observed, computed or inferred, and which endpoint a result actually measured.

- **Who (hypothesis):** families, students and new researchers; not known until asked.
- **Smallest artifact:** one Markdown page with each item tied to the rule it comes from and the ledger entry or reading-guide line that motivated it.
- **Check:** a different model reads the page beside the guardrails and flags any item the source does not support; each flag is confirmed on the source by a person.
- **Out of scope:** telling anyone what a paper means for their care; ranking journals or treatments.
- **Infrastructure and first week:** none; draft, then a blind read.
- **Stop rule:** if a reader with CF or a carer says it reads as instruction about their care, rewrite or withdraw it.

### 3. Replication cards

Take one headline number from a public table or file whose terms permit it, recompute it from the raw rows, and publish the card: source, hash, the arithmetic written out, the result, and the difference if any. The model's own count of placed residues across four structure files is the shape of it.

- **Who (hypothesis):** researchers and students who want to see a number checked rather than quoted.
- **Smallest artifact:** one card and the few lines of standard-library code that produced it.
- **Check:** the published total in the source is the frozen reference; the card passes `check_numbers.py` against the source file; a second person reruns the code.
- **Out of scope:** survival figures, anything cross-registry, anything the publisher does not allow us to reproduce.
- **Infrastructure and first week:** a script per card; one card from a structure file, where we already know the data.
- **Stop rule:** if a card needs a quotation over the overlap limit or a source marked manual or forbidden, it is not a card.

### 4. A negative-results shelf

A short page of things the lab tried and stopped, with the count that stopped them: six attempts at a WebGL page from one packet, a relay pilot that added nothing, a second local model that found no new defects.

- **Who (hypothesis):** other small labs deciding where not to spend time.
- **Smallest artifact:** a table with one row per stopped idea: what, when, the measure, where the record is.
- **Check:** every count on the page resolves to a line in a skill or proposal; `check_numbers.py` against those files.
- **Out of scope:** blaming a model or a vendor; general claims about what models can do.
- **Infrastructure and first week:** none; three rows from the records named above.
- **Stop rule:** if a row cannot be traced to a saved record, it leaves the page.

### 5. An open question board

A list of questions people ask about CF research, each mapped to the public sources that could answer it and the access class of each source, from the catalog. No answers on the board, only routes.

- **Who (hypothesis):** community writers and newcomers, and the lab when it decides what to read next.
- **Smallest artifact:** a Markdown table: question, sources, access class, what the answer would still not establish.
- **Check:** every source named is a catalog id; a script fails the page if one is not.
- **Out of scope:** questions about an individual's genotype, dose or prognosis; those are routed to a care team, not a source.
- **Infrastructure and first week:** a script that checks ids; ten questions from our own notes.
- **Stop rule:** if questions arrive that are about a person, the board needs a gate before it grows.

### 6. A record behind every number we already published

Walk the lab's own public notes and attach an evidence record or a file line to each number, and mark or remove the ones that have neither.

- **Who (hypothesis):** mainly us, and anyone checking us.
- **Smallest artifact:** a per-file list of numbers with their source lines, generated by the number checker.
- **Check:** `check_numbers.py` passes on every note against narrow evidence; the overlap checker passes.
- **Out of scope:** rewriting the notes' arguments.
- **Infrastructure and first week:** the two checks that exist; one note per evening.
- **Stop rule:** none needed; it is a finite audit.

### 7. Wider evidence coverage: EU records

Add EU regulatory records and the EU trials register to the evidence loop, so an approval or a trial can be checked in more than one jurisdiction.

- **Who (hypothesis):** European researchers and educators; the loop currently reaches US approval records only.
- **Smallest artifact:** one new provider following the skill's own five-step recipe, after a person reads the publisher's terms and records them in the catalog.
- **Check:** a saved fixture and a parse test, like the eight providers already there; the conduct tests still pass.
- **Out of scope:** any source whose terms have not been read; scraping a page.
- **Infrastructure and first week:** a script and a catalog entry; read the terms and write the limitation sentence first.
- **Stop rule:** terms that forbid automated access end the idea for that source; link to it by hand instead.

### 8. What changed between registry report years

A note per registry listing which definitions, equations and sections changed between the 2020 and 2024 reports, built from the reading guides, each line with its report, year and page.

- **Who (hypothesis):** researchers comparing figures across years; not known whether a registry already publishes this.
- **Smallest artifact:** one note for one registry.
- **Check:** every quotation is a byte-exact substring of the report, checked by script; the overlap check passes at the registry's limit; a registry scientist reads it before it leaves draft.
- **Out of scope:** any comparison of values across years or registries; survival methodology.
- **Infrastructure and first week:** none beyond the checks; the reports are downloaded by hand where the catalog says so.
- **Stop rule:** a registry asks us to stop, or the overlap check fails and cannot be fixed by paraphrase.

## 3. Sparks that need little or no infrastructure

Shorter ideas, imaginative but built only from what the lab can already do: evidence records, the registry reading guides, the CFTR structure page and its data, the ASCII and SVG skills, a local worker with a verifier, and plain files in git. Each says what could go wrong and what it asks of nobody: no patient data, no one's story.

**The limitation line as the headline.** The spark: every evidence record ends with what it does not establish, and that line is usually the most honest sentence in the room. Print only those lines from a run's ledger as a short reading list on uncertainty. Weekend version: a script that reads a JSONL ledger and writes a Markdown list of `limitations` with the question each came from. What could go wrong: a reader takes a limitation as a verdict against a paper. What it asks of nobody: the records describe documents, not people.

**One residue, one still.** The spark: the model page already knows, for each of the 1,480 positions, what UniProt annotates and which of the four files places it. Draw a single still SVG for one residue a week, with a caption generated from the data that says it is one arrangement, under one set of experimental conditions, and nothing more. Weekend version: reuse `make_page_blocks.py` and the SVG skill for one residue, with reduced motion respected because nothing moves. What could go wrong: a variant position is read as meaning something for a person; so show topology and placement only, never a classification. What it asks of nobody: the data is a public protein file.

**The chain strip in text.** The spark: the page shows which positions each structure places as a strip; a screen reader cannot see it. Render the same strip as a line of text, dots for placed and gaps for absent, with residue 508 named, and run it through the ASCII skill's accessibility checker. Weekend version: one script, four lines of output, one caption. What could go wrong: a line that wraps on a phone and misreads; test it at phone width. What it asks of nobody: nothing; it is the same data in another form.

**The DOI roll call.** The spark: the repo cites papers; papers get corrected and retracted. Before each release, run `retractions` over every DOI in the public notes and commit the records. Weekend version: a grep for DOIs piped to one command with `--ledger`. What could go wrong: treating `clean=true` as proof; the record itself says absence of a notice is strong but not conclusive. What it asks of nobody: the DOIs are already public.

**A not-found diary.** The spark: `not_found` and `out_of_scope` are results, and a list of the questions that got them is a map of what the public record does not say, or what we did not know how to ask. Weekend version: a Markdown page appended by hand after each session, with the command run and the status. What could go wrong: reading "not found" as "does not exist"; the page says so at the top. What it asks of nobody: our own queries only.

**Two readers, one draft.** The spark: a person and the local worker each list, independently, which sentences of a draft its quotation file supports. Publish only the disagreements, with what the source turned out to say. Weekend version: one note we have already written, the thinking profile, a draft piece kept to about 10,000 characters, and a table of disagreements. What could go wrong: the model's flags are treated as findings; each is a lead confirmed on the source by a person. What it asks of nobody: the draft is ours.

**A permission map.** The spark: the catalog's five access classes are a lesson in what "open" means, and 53 rows are hard to feel. Draw them as a small grid, generated from the catalog, one glyph per source, grouped by class, with the counts beneath. Weekend version: a script using the ASCII canvas library; regenerate it whenever the catalog changes and fail the build if the counts in the caption drift from the file. What could go wrong: a stale picture; so the caption is generated, never typed. What it asks of nobody: the catalog is public already.

**What a still cannot show.** The spark: the model page lists what it leaves out, no pore, no gate, no medicine, no motion, and each absence is a small lesson in how structures are made. Write one card per absence: why a still picture of a protein cannot show that thing, with a pointer to a primary source. Weekend version: four cards in one Markdown file, read blind by a second model against the structure note. What could go wrong: a card drifts into describing how a medicine works; the rule is that the cards explain the limit of a picture, never a treatment. What it asks of nobody: the subject is a method, not a person.

**The ledger, re-tallied.** The spark: the first ledger has 98 entries in one file; the interesting number is whether class C, the inferences beyond the evidence, falls as guards accumulate. Re-run `tally.py` after each session and keep the class counts as a small table over time, drawn as a text sparkline. Weekend version: a second column in the tally output and a dated table. What could go wrong: mistaking fewer entries for fewer mistakes when it means less looking; the table says which detectors ran. What it asks of nobody: the mistakes are ours.

## 4. What would make us stop or change course

- A publisher or registry asks us to stop using its text, or its catalog entry becomes `forbidden`: the idea that depends on it stops that day.
- A person with CF or a carer says a sentence lands as instruction about their care or uses them as content, and one edit cannot fix it: the artifact is withdrawn until it can.
- Review effort exceeds the useful work, measured in hours against the board's estimates, which today have no measured return on any of the 97 tasks: rescope or stop, and write the row for the negative-results shelf.
- An idea turns out to need patient data or an individual's genotype to be useful: it ends there.
- An existing tool already does it well: link to the tool and stop.

## 5. How to propose an idea

Open an issue or a pull request in this repository and quote the idea you are answering, or add a row here. A proposal must name who it could help, marked as a hypothesis until someone from that group says so; the smallest artifact; the check that could fail it; what it must never do; and a stop rule. Every number in it comes from a file or a command's output, with the source beside it.
