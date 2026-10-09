# A brief and a pilot protocol for a CF research leader

Status: proposed, written 2026-10-09. A draft for the maintainer to edit or drop. Nothing has been sent and no research leader has been contacted. Every factual statement about the lab is traced in section E to a file or a command; where one could not be traced it says so. Statuses in section B were checked against `origin/main` of each repository on 2026-10-09. Thresholds in section C are proposals to agree with a partner, not results.

## A. The brief

### The 30-second version

> We are a small, independent, open-source lab, built in donated time. We cannot cure anything, and we never read a person's genotype. What we can do is make gathering and checking public evidence cheaper. Our tools ask official public sources, print every value beside the field it came from, and fail when a claim has no source behind it. Give us a list of variant names you care about, and let your own team judge whether our sourced profiles save you time. People decide. The tools show their sources.

### The one-page leave-behind

> **Sensei Ewok's CF Lab: what we can and cannot do**
>
> We are a small, independent, open-source research lab, built in donated time. We are not affiliated with the Cystic Fibrosis Foundation, any university or any company. This is research, not medical advice.
>
> **What we cannot do.** We cannot cure anything. We do not interpret any person's genotype, recommend a medicine, or predict anyone's outcome. Clinical decisions belong with care teams.
>
> **What we can do.** Much research time goes into gathering public evidence and reconciling sources that disagree. Our tools aim to make that cheaper, and to show their working:
>
> - A variant profile prints what ClinVar states about one protein change, each value beside the API field it came from. When several records match, it prints them all and chooses none.
> - A brief puts several public sources side by side for one variant, drug or trial, and marks where a person should look.
> - A citation check confirms that every DOI, PubMed id and trial id in a note exists, and flags titles that do not match.
> - A claims checker fails any claim whose quote is not exactly in its source, or whose number is not in the quote.
>
> Our tools use only sources our public catalog says they may reach by an official API or a permitted fetch. Sources that ask for a person, or for written permission, are linked, not copied. CFTR2 is one of these.
>
> **What we propose.** A small pilot on your terms. You name the CFTR variants you care about, by public name only. We return a sourced profile for each. Your team scores each one against how you would do it. We agree the measures and the pass mark before anything runs, and we publish the result either way.
>
> **What we will never ask for.** Patient data of any kind. Your variant list stays yours unless you agree otherwise.
>
> **How to check us.** Our notes, tools, tests and mistake log are public in our research repository. The ledger lists the errors caught in our own sessions and what caught each one.

Word count of the leave-behind: 346, counted by script (section E).

## B. Five advanced scenarios, verified

Status words: **built** means the file exists at `origin/main` and a check passed when run on 2026-10-09; **prototype** means a partial working piece exists; **proposal** means nothing that does it exists yet. "Repo" names the repository; paths are relative to it.

| Scenario | Status | What exists, with its proof | What is not built |
| --- | --- | --- | --- |
| Sourced variant profiles | built | cf-research `tools/variant_profile/variant_profile.py`; its README says it prints what ClinVar states about one protein variant in one gene. Observed: `python -m unittest discover -s tools/variant_profile` ran 50 tests, OK, offline. cf-skills `cf-evidence-loop` `brief variant` (in `scripts/evidence/briefs.py`); SKILL.md lists `brief variant "<ClinVar query>"`. Its tests were not run here (no pytest on this machine); `tests/test_briefs.py` holds 16 test functions, counted, not run | Any source beyond ClinVar and PubMed; CFTR2 content (link only) |
| A living map of therapy trials | proposal | The parts exist: `brief trial <NCTId>` and `trials "<condition>"` in cf-skills `cf-evidence-loop` (SKILL.md section 1). A search of all four repositories for monitoring or a trial map found none | Scheduled monitoring, change detection, any map or page |
| A registry comparability ledger | proposal | The reading guides exist in cf-research only: `.claude/skills/cffpr-report-reading`, `ecfspr-report-reading`, `cf-genetics-reading`. None is in cf-skills at `origin/main` (searched the tree for registry, cffpr, ecfspr, genetics). `cf-research-context` says cross-registry arithmetic is out of scope until a registry scientist reviews it | The ledger itself; board task T-0011 (external review by a registry scientist) is still proposed |
| Structure view for variant groups | proposal | The site's CFTR page (`.ai` repo, `src/cf/cftr/index.html`, a draft) places one variant by position ("Place one variant by its position") and marks UniProt's annotated variant positions, "209 entries". Read with `git show origin/main:src/cf/cftr/index.html` | Grouping variants, or drawing a partner's list, on the structure |
| An evidence audit of "cure" claims | proposal | Its tools are built: cf-research `tools/claims/check_claims.py` (observed: 36 tests, OK, run from `tools/claims`) and cf-skills `cite-check` (`scripts/evidence/cite_check.py`; `tests/test_cite_check.py` holds 25 test functions, counted, not run). A search of all four repositories found no audit of cure claims | The audit itself: no set of cure claims has been gathered or checked |

Also observed on 2026-10-09: the source tools' tests (`tools/sources`, `python -m unittest`) ran 48 tests, OK.

### Observed in one run on one variant

`python variant_profile.py CFTR F508del --json`, run once from `tools/variant_profile`, 2026-10-09 (UTC), after reading `security-browsing`. Observed in one run on one variant; not a benchmark.

| Measure | Observed |
| --- | --- |
| Wall-clock seconds, start to exit | 3.63 |
| Exit code | 0 |
| Sources used | 1, ClinVar through NCBI E-utilities (catalog ids `ncbi-eutils` and `clinvar`, both `access: api`) |
| Records matched | 3, so the tool printed `ambiguous (3 records match)` and chose none |
| Claims checked | 0. The profile tool checks no claims; the claims checker is a separate step that was not run on this output |
| Requests | not printed. Our inference from the README: three (one esearch that found records, one esummary, one efetch) |

The output was kept in a scratch folder, not committed and not pasted here. No key was used and no contact address was sent.

## C. Pilot protocol (to agree with the partner before any run)

**Question.** For variants a research team cares about, does a lab profile save that team time without adding errors?

**Input.** The partner names 20 CFTR variants by public protein-level name (for example F508del). No patient data, no counts of people, no clinic names. The list stays the partner's; the lab publishes it only if the partner agrees in writing.

**What the lab returns, per variant, within a time agreed in advance.**

1. The `variant_profile.py` output (ClinVar via E-utilities).
2. The `brief variant` output (ClinVar and PubMed via E-utilities).
3. A short sourced note, a claims file for it, and the claims checker's output.
4. The `cite-check` output on that note.

**Sources.** Only sources the catalog marks `fetch` or `api`. For this pilot that is `ncbi-eutils` and `clinvar` (both `api`), plus `crossref` and `clinicaltrials-gov` (both `api`) for the citation check. CFTR2 is `manual` in the catalog: link only, nothing fetched or copied, until written permission is recorded there. Every request runs under the tools' own pacing and request caps.

**Metrics, defined before any run.**

| Metric | How it is measured |
| --- | --- |
| Claim error rate | Claims failing `check_claims.py`, plus claims the partner's reader marks unsupported by their quote, over all claims |
| Citation resolve rate | Identifiers `cite-check` marks `PASS`, over all identifiers. `UNRESOLVED` is reported apart from `NOT FOUND` |
| What was missed | Items in the partner's own curation of each variant that the profile lacks, as the partner lists them |
| Human minutes | Partner minutes to review a profile, against the partner's estimate of minutes to do the same work their usual way. Lab minutes reported too |

**Proposed thresholds, to agree with the partner.** Success: claim error rate of 0 after the checker, and at most 1 in 20 by the partner's reading; every cited identifier resolves or is removed; no item the partner calls essential is missed; partner review time below their usual time for most variants. Failure: any of these not met. These are proposals, not results.

**Stop rules.**

- Any patient data appears anywhere: stop, delete it, tell the partner.
- A source answers `blocked` or rate-limited: stop using it for the pilot; never route around it.
- The first five profiles exceed the agreed error threshold: stop and review with the partner before going on.
- The partner asks to stop, at any time, for any reason.

**What is shared.** The profiles, the claims files and checker output, the scoring sheet, and a short results note written with the partner. **Not shared:** any patient data, ever; the partner's variant list, unless they agree; the partner's internal curation.

**What would count as the lab being wrong.** A claim the partner shows is not supported by its quote. A value that differs from ClinVar at check time. A citation that does not resolve. An essential item missed. Net human time lost. Each is published in the results note and in the lab's ledger.

## D. What not to say, and questions to expect

**What not to say.**

- No cure claims, and nothing that implies a timeline or a clinical benefit.
- No genotype interpretation, and no eligibility for any medicine.
- No implied affiliation or endorsement by any organisation.
- No "AI discovers". The tools copy fields from public records; people judge them.

**Questions to expect, with honest answers.**

| Question | Answer |
| --- | --- |
| Is this better than CFTR2? | No. CFTR2 is the CF-specific authority. We link to it and do not copy it |
| Does a language model write the facts? | The evidence tool calls no language model. Notes may be drafted with a model; each claim is then checked by script, and the partner reads each one |
| Has anyone outside the lab reviewed it? | Not yet. That is what the pilot is for |
| What happens when ClinVar submitters disagree? | The profile prints every submitter's classification and a "submitters disagree" flag. It never resolves the disagreement |
| Who pays, and who owns the result? | We propose no fee. The tools are open source. The partner decides whether their list and scores are published |
| What if your tool is wrong? | We say so in the results note and the ledger, with what caught it |

## E. Statements about the lab and what shows them

| Statement | File or command that shows it |
| --- | --- |
| Small, independent, open-source lab, built in donated time | cf-lab `.claude/skills/cf-research-context/SKILL.md` ("a donated-time, open-source research effort"); cf-research `LICENSE`, `tools/LICENSE` |
| Not affiliated with the Cystic Fibrosis Foundation | cf-research `proposals/2026-10-04-youtube-about.md` and `proposals/2026-10-08-registry-permission-letters.md` (the lab's stated wording) |
| Not affiliated with any university or company | Not traced to a file. A statement of position for the maintainer to confirm before use |
| Research, not medical advice; no genotype interpretation, no medicine recommended | cf-research `README.md` (it says the lab does not interpret anyone's genotype or recommend or dose any medicine) |
| The profile prints each value beside its API field, prints every match and chooses none | `tools/variant_profile/README.md`; observed in the one run above (`ambiguous (3 records match)`) |
| A brief puts sources side by side and marks where a person should look | cf-skills `.claude/skills/cf-evidence-loop/SKILL.md`, Briefs (`CHECK`: a person should look) |
| The citation check confirms DOIs, PubMed ids and trial ids exist and flags titles | Same SKILL.md, Cite-check |
| The claims checker fails a quote not in its source, or a number not in the quote | `tools/claims/README.md`, rules R2 and R3; 36 tests OK, run 2026-10-09 |
| Tools reach only catalog sources marked `api` or `fetch` | cf-skills SKILL.md section 3 ("Catalog is the only permission"); `tools/variant_profile/README.md`, Limits |
| CFTR2 is linked, not copied | `sources/catalog.yaml` entry `cftr2`: `access: manual`, "link-only" |
| The evidence tool calls no language model | cf-skills SKILL.md section 1 ("never calls a language model") |
| Notes, tools, tests and the mistake log are public; the ledger lists errors and what caught each | cf-research `README.md` (its ledger paragraph: mistakes caught in working sessions, and "what caught each one") |
| No outside review yet | cf-research `README.md` ("No registry scientist has reviewed the guides yet"); cf-lab `tasks/BOARD.md` T-0011 still proposed |
| Leave-behind length | Counted with a short `python` split on whitespace over the quoted block, Markdown marks removed: 346 words, 2026-10-09 |

## Verification

- Statuses: `git show` and `git ls-tree` on `origin/main` of cf-skills, cf-lab and `.ai`; this worktree is cf-research `origin/main`.
- Tests run 2026-10-09: `tools/variant_profile` 50 OK; `tools/claims` 36 OK; `tools/sources` 48 OK. cf-evidence-loop tests not run: pytest is not installed here, and installing it was out of scope.
- Checks on this file: `tools/sources/check_numbers.py` with `proposals/2026-10-09-research-leader-brief-numbers-allow.txt`, and `tools/sources/check_source_overlap.py --quotes-only`.
