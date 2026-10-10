# Small AI tools for CF research

**Proposal inventory, refreshed October 6, 2026.** The first version (October 4, 2026) named three starting points and one experiment. This version says what each became, what we measured, scores what we could do next, and lists sparks that need almost no machinery. These are hypotheses about useful work, not validated community needs or clinical products. Check existing alternatives and seek qualified feedback before choosing a project. Research, not medical advice.

Every count below was read from a file or a command's output on 2026-10-06; the table listing each with its source is kept privately.

## 1. Where we are

### What the October 4 proposals became

| Proposal | Status | Evidence (observed) |
| --- | --- | --- |
| Literature evidence cards | Partly shipped, in a different shape | `cf-evidence-loop` (skills repo): 15 commands, 14 reaching one of 7 public hosts and one listing the catalog. A record carries 10 fields in the code (7 documented in the skill's reading table) and is refused with an empty `limitations`. Not shipped: the study-level card and the frozen-reference comparison; the 2026-10-05 local-worker evaluation records that no paired run against a cloud model was made |
| Historical claim ledger | Not started as proposed | No sourced timeline exists in the repo. Instead: an error ledger of our own mistakes (98 entries, `ledger/tally.py`) and three reading guides under `.claude/skills/`, two recording how a registry's definitions and layout changed between its 2020 and 2024 reports, one how a variant list's classes changed |
| Research-tool catalog | Partly shipped | `sources/catalog.yaml`: 53 sources and what each publisher permits (15 fetchable after `robots.txt`, 14 APIs with a recorded policy, 20 by hand, 3 refused, 1 by written request). `landscape/cf-repositories.md`: 6 repositories with activity and licence metadata. No comparison of data-access terms |
| First experiment | Not run as designed | What ran were loop measurements on code tasks (below), not frozen-field evidence cards |

### What we learned (measured; small samples; one lab, one machine)

1. **Our own re-read caught most mistakes; the costly class was caught by diversity.** Of 98 ledger entries, 86 have a classified first detector and 12 do not. Among the 86: our own re-read 45, an interpreter or test 20, a script assertion 12, a review by a frontier model 7, a person 2. Class C (inference beyond the evidence) was the largest, 43 entries, caught mainly, the ledger says, by a different model reading blind.
2. **A quiet review is not coverage.** On an 808-line script with 17 confirmed defects, four bare-packet runs of the local worker found 1, 0, 2 and 2; across both arms (8 runs), 9 of the 17 were found by no run.3. **Verified steps worked where one shot did not.** A page that failed 0 of 6 attempts in one packet passed 1 of 3 runs retried whole and 3 of 3 chains as five verified steps (arms unmatched on calls, 7 against 3; one task).
4. **Mode and data boundary each changed what the worker followed.** Fast mode followed 9 of 16 planted instructions in a plain prompt and 0 of 16 with the text wrapped as untrusted data; thinking mode followed 0 of 16 without any boundary.
5. **Rephrasing a request did not help.** A four-task, three-pass relay pilot: direct local work 3 of 12 accepted, one preparer plus local 3 of 12, another 0 of 12, over 60 local calls.

### What is still hypothesis

Who would use any of this, and whether an existing tool already serves them: not known; to be asked of the community. No registry scientist has read the reading guides. The CFTR model page stays a draft until an outside scientific review and a CF-community read (board task T-0081, proposed). The board holds 97 tasks, 96 proposed by a model, 1 by a person, none with a measured return.

## 2. Scored ideas

**Scoring (inferred; one author's judgement):** four marks from 1 to 5, summed: reach (kinds of reader), trust (how much a stranger can check us), effort (5 a weekend, 1 weeks), risk (5 almost no way to mislead, 1 a real chance). A sum is not impact or value per effort; it is four opinions added. Equal sums are tied, in reading order.

| Place | Idea | Reach | Trust | Effort | Risk | Sum |
| --- | --- | --- | --- | --- | --- | --- |
| 1= | How to read a CF paper, a checklist | 4 | 3 | 5 | 4 | 16 |
| 1= | Replication cards | 3 | 5 | 4 | 4 | 16 |
| 1= | A negative-results shelf | 2 | 4 | 5 | 5 | 16 |
| 4= | Source lookup before publishing | 4 | 5 | 3 | 3 | 15 |
| 4= | An open question board | 3 | 3 | 4 | 5 | 15 |
| 4= | A record behind every number we already published | 2 | 5 | 3 | 5 | 15 |
| 7 | Wider evidence coverage: EU records | 3 | 4 | 2 | 3 | 12 |
| 8 | What changed between registry report years | 3 | 4 | 2 | 2 | 11 |

Source lookup gets reach 4, not 5: its output is a raw evidence record, the most technical artifact here. By audience (hypotheses): a family member, the checklist then the question board; a community writer, source lookup then the checklist; a CF researcher, replication cards, registry report years and EU records (low here only on effort). A second model's agreement is never a reference answer.

**Check what exists first** (named, not endorsed): for the checklist, the Critical Appraisal Skills Programme (CASP) checklists and Testing Treatments, a public book on why treatments need testing; for replication cards, CODECHECK, which certifies independent execution and says its checkers record but do not investigate or fix; for source lookup, the Retraction Watch Database, a Crossref partner the tool's `retractions` command already reads; for the question board, the James Lind Alliance Priority Setting Partnership in CF (Rowbotham and colleagues, Thorax 2018, volume 73, issue 4, pages 388 to 390), which set research priorities with people with CF and clinicians, where ours would map questions to sources. If ours adds nothing to these, link to them and stop.

### 1= How to read a CF paper, a checklist

One plain-language page from our guardrails: the denominator, the year on the title page, observed or inferred, which endpoint a result measured.

- **Who (hypothesis):** families, students, new researchers.
- **Artifact and check:** each item tied to its rule and ledger line; a different model flags unsupported items, each confirmed on the source by a person. No infrastructure.
- **Out of scope; stop rule:** what a paper means for anyone's care; ranking journals or treatments. If a reader with CF or a carer says it reads as instruction about their care, rewrite or withdraw; until such a reader is found, it stays draft.

### 1= Replication cards

Recompute one headline number from a public table or file whose terms permit it; publish source, hash, the arithmetic written out, and the result.

- **Who (hypothesis):** researchers and students.
- **Artifact and check:** one card and the standard-library lines behind it; the source's published total is the frozen reference; `check_numbers.py` against the source file; a second person reruns it. The first card would recount the model's own count of placed residues across four structure files, our derived number; a true replication needs a catalog source with public raw rows and permissive terms, none named yet.
- **Out of scope; stop rule:** survival figures; anything cross-registry. A card needing a quotation over the overlap limit, or a manual or forbidden source, is not a card.

### 1= A negative-results shelf

Things the lab tried and stopped, with the count that stopped them: the relay rephrasing pilot (3 of 12 against 3 of 12 and 0 of 12); a second local model as reviewer, which added no defects; a written review card for the worker, which showed no claimable benefit.
- **Who (hypothesis):** other small labs; null-results journals exist, this shelf is lab-internal.
- **Artifact and check:** one row per stopped idea (what, when, measure, record); every count resolves to a line in a skill or proposal via `check_numbers.py`. No infrastructure.
- **Out of scope; stop rule:** blaming a model or a vendor; general claims about models. An untraceable row leaves the page.

### 4= Source lookup before publishing

Paste one sentence you intend to publish; get back the evidence record for the source it names, or an honest "not found". **A record never certifies the sentence**: it does not say the sentence is correct, complete or safe health information, it is not clinical guidance, and an approval record is a date, not a statement about who is eligible.

- **Who (hypothesis):** community writers, educators, the lab.
- **Artifact and check:** a documented mapping from sentence patterns (a DOI, an NCT id, a listed brand) to one of the 15 commands, printing `limitations` first; no model turns free text into a command. A frozen set of 30 sentences with known records, written before any run; count records returned, wrong records, sentences correctly declined. It makes network calls, bounded by the skill's rules: catalog permission only, paced requests, at most 200 attempts a process, cooldown on refusal.
- **Out of scope; stop rule:** judging whether the sentence is true; sentences about a named person; anything the catalog does not permit. If 3 or more of the 30 frozen sentences return a confident-looking wrong record, stop and redesign the output so "not found" is unmistakable.

### Ideas 4= to 8, in brief

| Idea and who (hypothesis) | Artifact and check | Out of scope; stop rule |
| --- | --- | --- |
| 4= Open question board: questions mapped to the public sources that could answer them and each source's access class; routes, not answers. Writers, newcomers, us | Table: question, sources, access class, what the answer would still not establish; a script fails the page if a source is not a catalog id | No question about an individual's genotype, dose or prognosis. Stop if questions about a person arrive before a gate exists |
| 4= A record behind every number we already published. Us, and anyone checking us | Per-note list of numbers with source lines; `check_numbers.py` and the overlap checker pass on every note | A finite audit |
| 7 Wider coverage: EU regulatory records and the EU trials register. European researchers and educators; the loop reaches US approval records only | One provider by the skill's five-step recipe after a person reads the terms, with a fixture and parse test | No unread terms, no scraping. Stop when terms forbid automated access |
| 8 What changed between registry report years. Researchers comparing years; not known whether a registry publishes this | One note per registry from the reading guides, each line with report, year, page; byte-exact quotation and overlap checks; a registry scientist reads it first | No comparison of values across years or registries. Stop if a registry asks |

## 3. Sparks that need little or no infrastructure

Built from what exists: evidence records, the structure data, the ASCII and SVG skills, plain files in git. Each names a wrong turn; none asks for patient data or anyone's story.

**The limitation line as the headline.** Every record ends with what it does not establish. Weekend: list the `limitations` lines from a run's ledger with their questions. Wrong turn: a limitation read as a verdict against a paper.

**One residue, one still.** For each of the 1,480 positions the page knows what UniProt annotates (CC BY 4.0, credited) and which of the four files places it. Weekend: one still SVG for one residue via `make_page_blocks.py`, captioned from the data as one arrangement under one set of conditions. Wrong turn: a variant position read as meaning something for a person; show topology and placement only.

**The chain strip in text.** Which positions each structure places, as dots and gaps a screen reader can read, residue 508 named. Weekend: one script, one caption, checked with the ASCII skill's accessibility checker from the control repo. Wrong turn: wrapping on a phone.

**The DOI roll call.** Before each release, run `retractions` over every DOI in the public notes and commit the records, within the 200-request budget. Wrong turn: `clean=true` read as proof; absence of a notice is strong but not conclusive.

**A permission map.** The 53 sources as a small grid by access class, generated from the catalog with the ASCII canvas library (control repo); a check run by hand before commit fails if the caption's counts drift from the file, since this repo has no CI. Wrong turn: a stale picture.

**What a still cannot show.** One card per absence the model page declares (no pore, no gate, no medicine, no motion): why a still picture cannot show it, with a primary source; read blind by a second model against the structure note. Wrong turn: drifting into how a medicine works; the cards explain the limit of a picture, never a treatment.

## 4. Boundaries

- Public aggregate information or synthetic fixtures only. No identifiable health data in repositories, logs or prompts.
- Contributors and reviewers are never asked whether they or someone close to them has CF, nor to disclose personal medical information.
- Notes stay toolchain-free; code lives under `tools/`, one directory per tool, with its own tests.
- Treat [CFF priorities](cff-goals.md) and the [existing-tool landscape](landscape/medical-ai-landscape.md) as starting evidence, not proof that any of this is needed.

## 5. What would make us stop or change course

- A publisher or registry asks us to stop using its text, or its catalog entry becomes `forbidden`: the dependent idea stops that day.
- A reader says a sentence lands as instruction about their care or uses people as content, and one edit cannot fix it: withdrawn until it can.
- Review effort exceeds the useful work, measured in hours (no board task has a measured return yet): rescope or stop.
- An idea needs patient data or an individual's genotype to be useful: it ends there.
- An existing tool already does it well: link to it and stop.

## 6. How to propose an idea

Tell the lab on the Contact page of [senseiewok.ai](https://senseiewok.ai/contact/); this repository does not take outside pull requests or issues. A proposal names who it could help, marked as a hypothesis until public feedback from people in that role says otherwise, with nobody asked about their own health; the smallest artifact; the check that could fail it; what it must never do; a stop rule. Every number comes from a file or a command's output, with the source beside it.
