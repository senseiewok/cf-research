# External Permission Letters & Outbound Dispatch Dossier

Status: Prepared 2026-10-16 for maintainer dispatch (Task T-0106; gates Section 2 of `tasks/HUMAN-TASKS.md`).

This dossier unifies all external permission and review request letters into one actionable dispatch document. Per lab governance and `tasks/HUMAN-TASKS.md` Section 2, **all letters must be reviewed, filled, and sent by the human maintainer from the lab's official email address (`[lab contact address]`)**. AI agents prepare, standardize, and verify drafts, but never send outbound emails from the lab address or assert permissions without written replies.

---

## 1. Dispatch Checklist for Maintainer

Check off each letter as sent, noting the dispatch date and recipient address used:

- [ ] **Letter 1: CFTR2 Permission Letter** (Johns Hopkins University & CF Foundation)
- [ ] **Letter 2: Child Opportunity Index (COI) Letter** (diversitydatakids.org, Brandeis University)
- [ ] **Letter 3: CF Foundation Patient Registry Letter** (CFF Registry Team)
- [ ] **Letter 4: ECFSPR Annual Data Report Letter** (European Cystic Fibrosis Society)
- [ ] **Letter 5: UK Cystic Fibrosis Registry Letter** (Cystic Fibrosis Trust)
- [ ] **Letter 6: Australian CF Data Registry (ACFDR) Letter** (Monash University)
- [ ] **Letter 7: Clinical & Care-Team Review Invitation** (AI-in-care safe-use review)

---

## 2. Reply Tracking Ledger

When replies arrive (or if no response is received after 60 days), record the outcome here and open a companion pull request updating `sources/catalog.yaml`:

| Target Source | Catalog ID | Current Status | Date Sent | Reply Date | Outcome | Permitted Uses / Scope |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **CFTR2** | `cftr2` | Drafted, awaiting maintainer dispatch | Pending | Pending | Pending | Link-only until written consent |
| **Child Opportunity Index** | `coi-diversitydatakids` | Drafted, awaiting maintainer dispatch | Pending | Pending | Pending | Link-only until terms clarified |
| **CFF Patient Registry** | `cffpr-adr-2024` | Drafted, awaiting maintainer dispatch | Pending | Pending | Pending | Manual read, link-only |
| **ECFSPR** | `ecfspr-adr-2024` | Drafted, awaiting maintainer dispatch | Pending | Pending | Pending | Fetch allowed; redistribution pending |
| **UK CF Registry** | `ukcfr-adr-2024` | Drafted, awaiting maintainer dispatch | Pending | Pending | Pending | Short quotes & single values pending |
| **Australian CFDR** | `acfdr-adr-2025` | Drafted, awaiting maintainer dispatch | Pending | Pending | Pending | Single cited values pending |
| **AI-in-Care Review** | N/A (Clinical safety) | Drafted, awaiting maintainer dispatch | Pending | Pending | Pending | Independent clinical validation |

---

## 3. General Sending Rules

1. **Sender Address**: Use the official lab address (`[lab contact address]`), never a personal one or personal name; the lab speaks for itself as an open-source research initiative.
2. **Fill All Placeholders**: Fill in `[sender name]`, `[lab contact address]`, `[date]`, and recipient addresses.
3. **No Private Material**: Never include patient data, API keys, or private internal notes.
4. **Silence Is Not Permission**: If no response arrives, the lab strictly adheres to current catalog terms: read and cite, link-only, no re-hosting, no bulk redistribution.
5. **Private Reply Storage**: Store full reply emails in untracked `../cf-lab-files/` (never committed to git). Record public outcomes and dates in `sources/catalog.yaml` and this ledger.

---

## 4. Letter Templates

### Letter 1: CFTR2 (Variant Annotation & Classification)

- **Recipient**: Address specified on publisher's terms page (`info@cff.org` per legal terms PDF)
- **Catalog Entry**: `sources/catalog.yaml` (`cftr2`)
- **What We Ask**: Written permission to display variant name, clinical classification, and functional classification for single searched variants on `https://senseiewok.ai/cf/cftr/`, with unaltered wording and full attribution.

> **Subject**: Request for written permission to show a few CFTR2 fields on an educational page
>
> Hello,
>
> We are a small, independent, open-source project that runs on donated time. We are not affiliated with the Cystic Fibrosis Foundation, Johns Hopkins University, any conference or any company.
>
> We publish a draft, educational 3D model of the CFTR protein at https://senseiewok.ai/cf/cftr/. A visitor can type a protein-level variant name and see where that position sits on the structure; the page uses none of CFTR2's data and says only that it shows a position. Your terms say that republishing or distributing any portion of CFTR2 needs written permission, so we are asking first. Today we only link to CFTR2.
>
> We would like to show, for one variant a visitor types, these fields only:
> - the variant name;
> - the clinical classification, as CFTR2 states it;
> - the functional classification, as CFTR2 states it;
> - an attribution line and a link to that variant's CFTR2 page;
> - the date we retrieved the data.
>
> The data would come from a download or method you agree to, at a pace you set. CFTR2's wording would appear unaltered.
>
> What we would not do:
> - host or share anything beyond the agreed fields, in bulk or otherwise;
> - interpret anyone's genotype, or offer medical advice;
> - collect or show patient data;
> - use the data commercially;
> - change your wording.
>
> The page would say plainly that it is an independent, educational project, not affiliated with or endorsed by CFTR2, the Foundation or Johns Hopkins, and that visitors should talk with their care team.
>
> What we offer: credit in the wording you choose, a link back on every result, a report to you of any data error we notice, and open access to our source code and tests.
>
> If the answer is no, or if we do not hear back, that is fine. We will keep linking to CFTR2 and show none of its data.
>
> Thank you for building and maintaining CFTR2.
>
> [sender name], for the lab · [lab contact address] · [date]

---

### Letter 2: Child Opportunity Index (COI)

- **Recipient**: Address from diversitydatakids.org contact page
- **Catalog Entry**: `coi-diversitydatakids`
- **What We Ask**: Confirmation on using published aggregate neighbourhood scores in an open educational resources guide.

> **Subject**: Using Child Opportunity Index scores in an open educational guide
>
> Hello,
>
> We are a small, independent, open-source project that runs on donated time. We are not affiliated with Brandeis University, any foundation, any conference or any company.
>
> We are writing an open, free guide on how people and care teams can find social-needs resources, such as help with food, housing or transport. We would like to use the Child Opportunity Index's neighbourhood scores in it, and we would rather ask than guess.
>
> Two questions:
> 1. May we use the index's published neighbourhood scores in this guide, and are there conditions?
> 2. How would you like the index cited?
>
> What we would do: use aggregate neighbourhood scores only, cite the index and its version, link to your site, and explain in plain words what the scores measure.
>
> What we would not do: never place an individual's address, or any child, on the index; never collect personal data; never use the scores to judge a person or a family; no commercial use. The guide would say clearly that the index describes a neighbourhood, not a family.
>
> If the answer is no, or if we do not hear back, we will link to your site and not use the scores.
>
> Thank you for making this work public.
>
> [sender name], for the lab · [lab contact address] · [date]

---

### Letter 3: Cystic Fibrosis Foundation Patient Registry (CFFPR)

- **Recipient**: Registry permissions email from Annual Data Report Figure Permissions notice
- **Catalog Entry**: `cffpr-adr-2024`
- **What We Ask**: Permission to extract published aggregate numbers into a machine-readable table citing report, year, and page with stated denominators; programmatic access preferences.

> **Subject**: Question about using published Annual Data Report figures with citation
>
> Hello,
>
> We are a small, independent, unfunded, open-source research project. We build tools that help people check claims about cystic fibrosis against published sources. We are not affiliated with the Foundation and do not claim its endorsement.
>
> We would like to read the published Patient Registry Annual Data Reports and keep a table of the aggregate numbers they print, where every value carries the report title, data year, page and web address, together with the denominator the report states. No patient-level data is involved, and none is being requested.
>
> Three questions:
> 1. Is it acceptable to extract numeric values from the published reports and share them in such a table, with that citation?
> 2. What attribution wording would you like?
> 3. Automated requests to cff.org are refused, so we download the reports by hand. Is there an arrangement you prefer for programmatic access, or should we continue by hand?
>
> For transparency: software, including AI language models running on our own computers, helps read the reports. We send no report text to any outside service unless you tell us that is acceptable. If the answer to question 1 is no, we will publish only the extraction method and no values.
>
> Thank you for the registry and for publishing the reports openly.
>
> [sender name], for the lab · [lab contact address] · [date]

---

### Letter 4: European Cystic Fibrosis Society Patient Registry (ECFSPR)

- **Recipient**: ECFSPR registry contact (`ecfspr@ecfs.eu` per annual report front matter)
- **Catalog Entry**: `ecfspr-adr-2024`
- **What We Ask**: Confirmation on redistributing aggregate cited values in public research tools and attribution wording.

> **Subject**: Reuse of published ECFSPR Annual Report figures with citation
>
> Hello,
>
> We are a small, independent, unfunded, open-source research project (not affiliated with ECFS or the registry). The annual reports are openly published and your site permits automated retrieval, so we read them with software.
>
> We would like to keep a table of the aggregate numbers the reports print, where each value carries the report title, data year, page, web address and the denominator the report states. No patient-level data is involved.
>
> 1. May we share such a table publicly, with that citation?
> 2. What attribution wording would you like?
> 3. Are any figures or tables excluded from reuse?
>
> For transparency: AI language models running on our own computers help read the reports; no report text goes to an outside service unless you say that is acceptable.
>
> [sender name], for the lab · [lab contact address] · [date]

---

### Letter 5: UK Cystic Fibrosis Registry (Cystic Fibrosis Trust)

- **Recipient**: Registry contact per copyright notice in 2024 Annual Data Report
- **Catalog Entry**: `ukcfr-adr-2024`
- **What We Ask**: Written permission to quote short phrases and cite single published values with page references.

> **Subject**: Permission to cite and briefly quote the UK CF Registry Annual Data Report
>
> Hello,
>
> We are a small, independent, unfunded, open-source research project (not affiliated with the Trust). The 2024 Annual Data Report says that using or reproducing its content in publications needs the Trust's permission, so we are asking before we do either.
>
> We would like to cite single values with the report's title, page and web address, and to quote at most a short phrase where needed. We would not reproduce tables or figures and would not host the report.
>
> 1. Do you give permission for that use?
> 2. What attribution wording would you like?
>
> For transparency: AI language models running on our own computers help read the reports; no report text goes to an outside service unless you say that is acceptable.
>
> [sender name], for the lab · [lab contact address] · [date]

---

### Letter 6: Australian Cystic Fibrosis Data Registry (ACFDR)

- **Recipient**: Monash University registry coordinator per 2025 Annual Report enquiries notice
- **Catalog Entry**: `acfdr-adr-2025`
- **What We Ask**: Confirmation that single cited values and short quotations are acceptable.

> **Subject**: Enquiry about citing the Annual Report with short quotations
>
> Hello,
>
> The 2025 Annual Report sends enquiries about use or reproduction to the registry's contact, so we are writing before we do either. We are a small, independent, unfunded, open-source research project (not affiliated with the registry).
>
> We would like to cite single published values with the report number, page and web address, and to quote at most a short phrase. We would not reproduce tables or host the report.
>
> 1. Is that acceptable?
> 2. What attribution wording would you like?
>
> AI language models running on our own computers help read the reports; no report text goes to an outside service unless you say that is acceptable.
>
> [sender name], for the lab · [lab contact address] · [date]

---

### Letter 7: Invitation to Review AI-in-Care Materials (Clinical Safety)

- **Recipient**: CF clinician, social worker, or person with CF / family carer (obtained via direct relationship)
- **Context**: Gated by `tasks/HUMAN-TASKS.md` Section 3 (T-0105 safe-use review)

> **Subject**: Would you review our AI-in-care page? About 30 minutes
>
> Hello,
>
> We are a small, independent, open-source project that runs on donated time, not affiliated with the Cystic Fibrosis Foundation, any conference, any university or any company. We have written guidance on using AI chat tools alongside CF care, and we would value the view of someone who lives or works close to it: a CF clinician, a social worker, or a parent or carer.
>
> The draft page is here: https://senseiewok.ai/resources/ai-in-cf-care/ (and in our open repository).
>
> What we ask: about 30 minutes of your time to read the page and tell us:
> 1. Does anything sound like medical advice, or risk someone skipping their care team?
> 2. Are the crisis and emergency contact recommendations clear and appropriate?
> 3. Does the tone treat people with dignity?
>
> What we do not want: we will never ask for information about any patient or family.
>
> If you have time to look, we would be grateful. If not, no reply is needed at all.
>
> [sender name], for the lab · [lab contact address] · [date]
