# Registry permission letters: drafts to send (board T-0006, T-0072)

Status: proposed, written 2026-10-08. These are drafts, not sent messages. Only the maintainer decides whether to send and sends each one, from a lab-owned address. Nothing is sent, and no registry has been contacted. The facts about each registry come from `sources/catalog.yaml` (entry ids in brackets), read 2026-10-04 to 2026-10-05; they are the lab's notes, not the registries' current terms, so re-read each report's own copyright notice before sending.

This replaces the single draft in `landscape/registry-data-access.md` with one letter per registry, because the four registries ask for different things. The older draft stays until a reply is recorded there.

## What the lab is asking, and why it differs by registry

| Registry | What its report says about reuse (from the catalog) | What we ask |
| --- | --- | --- |
| CF Foundation Patient Registry [`cffpr-adr-*`] | A Figure Permissions paragraph tells readers to contact the Registry team by email to request use of its charts and data. Automated requests to cff.org are refused, so reports are downloaded by hand | Permission to extract published numbers into a table that cites report, year and page; an attribution string; whether a programmatic route exists |
| ECFSPR [`ecfspr-adr-*`] | robots.txt permits automated retrieval. Its reuse terms are unread | Whether extracted values may be redistributed with citation, and the attribution string |
| UK CF Registry (Cystic Fibrosis Trust) [`ukcfr-adr-2024`] | The copyright notice says using or reproducing content in publications needs the Trust's permission | Written permission to quote short passages and to cite single values with page references |
| Australian CF Data Registry [`acfdr-adr-2025`] | States no licence; its only reuse statement sends enquiries to the registry's Monash contact | Whether short quotations and single cited values are acceptable |

## Before sending: facts to fill in or check

- Sender line: the lab's name and a lab-owned address. Never a personal address or the maintainer's name; the lab speaks for itself.
- Recipient: copy each address from the report's own permissions or contact paragraph. The catalog does not hold them, and `info@cff.org` is named only for CFTR2, not for the registry.
- Re-read the report's copyright paragraph and note its date in the reply log below.
- Do not attach or paste any report text beyond a short citation.

## Letter 1: CF Foundation Patient Registry

> Subject: Question about using published Annual Data Report figures with citation
>
> Hello,
>
> We are a small, independent, unfunded, open-source research project. We build tools that help people check claims about cystic fibrosis against published sources. We are not affiliated with the Foundation and do not claim its endorsement.
>
> We would like to read the published Patient Registry Annual Data Reports and keep a table of the aggregate numbers they print, where every value carries the report title, data year, page and web address, together with the denominator the report states. No patient-level data is involved, and none is being requested.
>
> Three questions:
>
> 1. Is it acceptable to extract numeric values from the published reports and share them in such a table, with that citation?
> 2. What attribution wording would you like?
> 3. Automated requests to cff.org are refused, so we download the reports by hand. Is there an arrangement you prefer for programmatic access, or should we continue by hand?
>
> For transparency: software, including AI language models running on our own computers, helps read the reports. We send no report text to any outside service unless you tell us that is acceptable. If the answer to question 1 is no, we will publish only the extraction method and no values.
>
> Thank you for the registry and for publishing the reports openly.
>
> [Lab name] · [lab-owned address]

## Letter 2: ECFSPR

> Subject: Reuse of published ECFSPR Annual Report figures with citation
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
> [Lab name] · [lab-owned address]

## Letter 3: Cystic Fibrosis Trust (UK Registry)

> Subject: Permission to cite and briefly quote the UK CF Registry Annual Data Report
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
> [Lab name] · [lab-owned address]

## Letter 4: Australian CF Data Registry (Monash)

> Subject: Enquiry about citing the Annual Report with short quotations
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
> [Lab name] · [lab-owned address]

## After sending

Record, in `landscape/registry-data-access.md`, the date sent, the recipient route (not the address), the date and gist of any reply, and what it permits. A reply is evidence for one use only. Until a reply arrives, the catalog's current rules stand: read and cite, do not re-host, do not reproduce tables. Do not send report text to a cloud model for a source whose terms are unread (board T-0072, point 1).

Not covered: the Canadian registry's terms and the CF Foundation's current Figure Permissions wording, which were not re-read for this draft; whether a registry will reply; and any legal view (a lawyer's opinion is worth having if the lab ever needs certainty).
