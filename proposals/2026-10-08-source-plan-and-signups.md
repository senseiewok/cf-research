# Bringing sources together: a ranked plan and the sign-ups it needs

Status: proposed, written 2026-10-08. The maintainer decides what to pursue and does every sign-up; nothing here is admitted, and no account, key or licence has been created. The ranking is the controlling agent's judgment of impact on one job: answering a CF question from several checked sources at once. The access facts marked "catalogued" come from `sources/catalog.yaml`; everything marked "unread" is general knowledge that has not been checked against the source's own terms.

## Ranked by impact

| # | Source | What it adds | Access | Status |
| --- | --- | --- | --- | --- |
| 1 | The sources already approved, wired together: ClinicalTrials.gov, openFDA, ClinVar, PubMed, Europe PMC, Crossref with Retraction Watch, OpenAlex, NIH RePORTER | One question gives one cross-checked record (trial, label, variant, paper, retraction status) | No key, except OpenAlex (key required) and optional keys for NCBI and openFDA | Catalogued |
| 2 | Variant layer: ClinVar (`tools/variant_profile`, open PR), MyVariant, AlphaMissense v3, CFTR2 | Evidence profiles for rare variants, plus an independent prediction score | Open; CFTR2 is a manual download | ClinVar and MyVariant catalogued; AlphaMissense in the draft catalog PR |
| 3 | Registry reports: ECFSPR, UK, Australia, CFRI; the CFF Patient Registry | Population numbers and trends | Four open; the CFF registry needs a permission request (board T-0006) | Catalogued |
| 4 | Drug labels: openFDA, DailyMed, EMA | Approval and label facts without model guesses | No key | openFDA and EMA catalogued; DailyMed in the draft catalog PR |
| 5 | Claude connectors for ChEMBL, ClinicalTrials.gov, PubMed and bioRxiv | Ready-made tools in the agent | Sign-in with the maintainer's Claude account | Listed in the session, not signed in. Their results do not pass through the lab's own audit trail |
| 6 | Structures: RCSB and UniProt (in use), AlphaFold DB | Variant overlays on the CFTR viewer | Open | AlphaFold unread |
| 7 | Open Targets and ChEMBL APIs | Gene-disease evidence and compound activity for modulators | Open | Unread |
| 8 | Tissue expression: Human Protein Atlas or Expression Atlas | Whether CFTR is expressed in a tissue, since GTEx's robots.txt blocked our agent | Open downloads | Unread |
| 9 | Paper trail: Semantic Scholar, Unpaywall, Scite | Citation context and free full text (board T-0044, T-0047) | Form, email or paid | Semantic Scholar and Unpaywall in the draft catalog PR |
| 10 | Population denominators: Data Commons (API and MCP) | Context for registry numbers (board T-0048) | Key required | In admission review |
| 11 | BioMCP | A discovery tool beside the evidence loop (board T-0042) | Admission review | Proposed |
| 12 | Additional model families as reviewers (Gemini and others) | A different family for blind reviews (board T-0101) | Key; may be paid | Proposed |
| 13 | GitHub MCP | Development tooling only (board T-0050) | Repo-scoped token | Proposed |

Not pursued: patient records, scraped pages, and any source the catalog marks forbidden (PubPeer, ClinGen, LOVD).

## What only the maintainer can do (in the order that unlocks the most)

Use an address that belongs to the lab, not a personal one, for every sign-up. Put any key only in the git-ignored `.env` in `cf-lab`, never in a chat, a file in a repo, or a pull request. Name the variable, never the value.

| Step | What you do | Where | What I do after |
| --- | --- | --- | --- |
| A | Create an OpenAlex API key (the catalog says a key is required for our use) and read its rate-limit page | docs.openalex.org, the terms URL in the `openalex` entry | Add `OPENALEX_API_KEY` to `.env.example`, build the provider, test with one request |
| B | Optionally create an NCBI key through an NCBI account (the catalog says optional; it raises the request limit, exact limit unread) | The NCBI account page linked from the E-utilities pages | Read the limit, add `NCBI_API_KEY`, raise the cap only to what the page states |
| C | Optionally create an openFDA key | open.fda.gov authentication page (the `openfda-label` terms URL) | Same as B |
| D | Send the CFF Patient Registry permission request, and the ECFSPR request (board T-0006) | The registries' contact pages; I draft the letters, you send them | Extract only what they permit, with the denominators the reports print |
| E | Sign in to the Claude connectors you want (ChEMBL, ClinicalTrials.gov, PubMed, bioRxiv) | Claude settings, connectors | Use them for discovery; anything that ships is re-checked through the catalog route |
| F | Read and accept Semantic Scholar's API licence for the lab if you want it (it binds the organisation; indemnity; Washington State law), then request a key through the form on its API page | The licence and overview pages in the `semantic-scholar-api` entry | Admit the entry and build the provider |
| G | Create a Data Commons key (free per the catalog's note) if you want denominators | apikeys.datacommons.org | Finish the admission review (T-0048) |
| H | Decide whether the Gemini route is worth a key, and whether your licence covers API use (board T-0101) | Google's AI developer pages | Build the one-packet adapter |
| I | Paid items (Scite Pro, any model API beyond a free tier): decide whether the value justifies it. I have not read any price | The vendor's pricing pages | Trial one seat only, per T-0047 |

## What stays true whatever you sign up for

- Every source goes into `sources/catalog.yaml` first, with the terms read and dated, and no `base_url` until a human admits it. A key never turns a forbidden source into an allowed one.
- Politeness caps stay: at most one request a second unless the source's own page states more, a descriptive user agent, no personal identifier.
- A source adds a new external network call and sometimes a new dependency; every pull request says so.
- Paid or key-bearing sources are one packet at a time: nothing private, nothing about a person, goes into a call. Approval to use a cloud source does not authorise sending patient or private material.
- Not covered here: the terms of every row marked unread, the cost of any paid item, and whether any subscription you already hold covers API use.
