# Proposal: sources a CF researcher actually uses, and which to add to the catalog

Status: proposal. Written 2026-10-04 from domain knowledge, not from fetching; every entry below is `[unverified]` for access terms until someone reads the publisher's terms and records the result in `sources/catalog.yaml`. The point of this note is to name the gaps precisely, so the access-checking work is bounded.

## The professional loop, and where the catalog stands

A working CF researcher answers a question by walking seven layers. The catalog's coverage today:

| Layer | What the researcher asks | Covered now | Gap |
| --- | --- | --- | --- |
| 1. Literature | What has been published? | Crossref, OpenAlex, NCBI E-utilities, Europe PMC (unverified), bioRxiv | Systematic reviews as a distinct, higher-evidence tier |
| 2. Integrity | Has it been retracted, corrected, or disputed? | Nothing | **Largest gap for a claim-verification lab** |
| 3. Trials | What is being tested, on whom, with what endpoint? | ClinicalTrials.gov v2 | Non-US registries; CF-specific trial networks |
| 4. Regulatory | What is actually approved, where, for whom? | Nothing (the website cites one DailyMed label by hand) | Labels and assessment reports as structured sources |
| 5. Population | How many people, with what outcomes? | CFFPR, ECFSPR | Other national registries that publish open aggregate reports |
| 6. Guidelines | What does the field recommend, and on what evidence grade? | Nothing | Care standards and consensus statements |
| 7. The field's own conversation | What was presented, funded, or argued before it was published? | Nothing | Conference abstract supplements, funding databases, preprint servers beyond bioRxiv |

Variant biology (CFTR2, CFTR-France, LOVD) sits across layers 1 and 5 and is already catalogued.

## Candidate sources by layer

`Access` is my expectation, not a test result. Record the tested value in the catalog.

### Layer 2: integrity (add first)

| Source | What it gives | Expected access | Why it matters here |
| --- | --- | --- | --- |
| Retraction Watch Database (via Crossref, open since 2023) | Retraction and correction records matched to DOIs; Crossref `update-to` relations | `api` through api.crossref.org, already catalogued | Every DOI in the registry table and every cited paper gets a mechanical retraction check. This is the single cheapest credibility feature the lab can ship |
| PubPeer | Post-publication comments by DOI | `api` with key, terms to read | Signals dispute, not verdict. Surface, never adjudicate |
| Cochrane Cystic Fibrosis and Genetic Disorders Group reviews | GRADE-rated systematic reviews on CF interventions | `manual` for full text; abstracts via PubMed | The field's own evidence grades. A claim that contradicts a current Cochrane CF review needs a very good reason |

### Layer 4: regulatory

| Source | What it gives | Expected access | Note |
| --- | --- | --- | --- |
| DailyMed / openFDA drug label API | Current US prescribing information, structured | `api` (openFDA is on the sandbox allowlist) | Turns the website's hand-cited label into a checked field with a version date |
| FDA Drugs@FDA | Approval history, dates, supplements | `api` via openFDA | Approval year claims (T-0023) become one lookup |
| EMA EPARs and product information | EU approvals, assessment reports | `manual`; terms to read | EU counterpart; needed before any "approved in Europe" sentence |
| MHRA, Health Canada, TGA product registers | Other major jurisdictions | `manual` | Lower priority; add when a claim needs them |

### Layer 5: population, beyond US and EU

| Source | Expected access | Note |
| --- | --- | --- |
| UK Cystic Fibrosis Registry annual data report (Cystic Fibrosis Trust) | `manual` likely; PDF | Third large registry; UK is part of ECFSPR but publishes its own, more detailed report |
| Canadian CF Registry annual report (Cystic Fibrosis Canada) | `manual` likely | Comparability notes needed, as for all cross-registry work |
| Australian CF Data Registry annual report | `manual` likely | Same |
| Cystic Fibrosis Registry of Ireland annual report | `manual` likely | Small; useful as a comparability case |

Do not add individual-level registry datasets. The "no patient data" line is absolute; aggregate published reports only.

### Layer 6: guidelines

| Source | Expected access | Note |
| --- | --- | --- |
| ECFS standards of care papers (Journal of Cystic Fibrosis) | Publisher terms; abstracts via PubMed | Primary European guideline series |
| CFF clinical care guidelines | `manual`; www.cff.org is 403 to clients | Same access problem as the registry reports |
| Cystic Fibrosis Trust consensus documents | `manual` likely | UK |

### Layer 7: the field's own conversation

| Source | Expected access | Note |
| --- | --- | --- |
| NACFC abstracts (published as a Pediatric Pulmonology supplement) | Publisher terms; indexed via Crossref | Observed 2026-10-04: Crossref already returns JCF supplement abstracts (e.g. `P226`, `EPS1.10` numbered entries), so conference material is reachable through an API we already use `[verified]` |
| ECFS conference abstracts (Journal of Cystic Fibrosis supplement) | Same | Same |
| medRxiv | `api` via api.biorxiv.org, already catalogued | Add explicitly; clinical preprints land here, not on bioRxiv |
| NIH RePORTER | `api`, public | Who is funded to do what; the best early signal of where the field is going |
| CFF-funded awards | `manual`; cff.org | Same access issue |
| EU Clinical Trials Register / CTIS | `api` or `manual`; terms to read | Non-US trials |
| ECFS Clinical Trials Network, CFF Therapeutics Development Network | `manual` | Context on how CF trials are run; not data sources |

### Variant and structural biology (already partly covered)

| Source | Expected access | Note |
| --- | --- | --- |
| ClinVar | `api` via E-utilities, already allowlisted | Variant clinical significance with submitter evidence; complements CFTR2 |
| gnomAD | `api` (GraphQL), terms to read | Aggregate allele frequencies; aggregate only, no individual data |
| UniProt P13569, RCSB PDB, AlphaFold DB | `api`; blocked by sandbox allowlist today, not by publishers | Needed only for the structural atlas, project #4 |
| CFTR1 (SickKids CF Mutation Database) | `manual`; terms to read | Historical nucleotide-level catalogue |

### Ontologies and identifiers (plumbing, not sources)

MeSH, MONDO, HPO, Orphanet. These make the catalog's indicator slugs and literature queries interoperable with everyone else's. Add when the extraction schema is frozen (T-0008), not before.

## What not to add

- Patient forums, Facebook groups, Reddit, X. `security-browsing` section 6 is explicit: themes, not people. These are never data sources.
- Individual-level datasets of any kind, including de-identified registry extracts and patient-derived omics (GEO, dbGaP). Out of scope by charter, regardless of what the data use agreement allows.
- News sites and press releases as evidence. Useful as leads to a primary source; never cited as one.

## Recommended order

1. **Retraction Watch via Crossref** — already permitted, already reachable, highest credibility per hour. Add to the catalog and to the claim-check harness design: every DOI the lab cites gets a retraction lookup, and the result is a field, not a footnote.
2. **openFDA label and approval APIs** — allowlisted, structured, resolves T-0023 and every future "approved in" sentence.
3. **medRxiv and ClinVar** — zero new permissions, both reachable through APIs already in the catalog.
4. **Cochrane CF reviews** — abstracts via PubMed now; full-text terms to read.
5. **NIH RePORTER** — public API, no key.
6. **UK, Canadian, Australian registry reports** — one afternoon of reading terms and recording access values; then the fetcher handles them like ECFSPR or like CFF.
7. Everything else as a claim needs it.

Each item is one catalog entry plus a terms check. None requires new code beyond the fetcher that exists.
