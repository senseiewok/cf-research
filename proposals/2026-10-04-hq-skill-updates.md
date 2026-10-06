# Proposed updates to HQ skills from the registry-access research

Status: **applied in hq on 2026-10-04 (local commit `f1fe8ff`, board task T-0003; not yet pushed or merged). It read "proposal, not applied" when written; the text below is the original proposal.** Files under `hq/.claude/` are write-protected from automated agents in this workspace, which is the right default. Each change below is the exact text to paste, with the evidence it rests on, so a human can apply it after reading. Per `AGENTS.md`: narrowly scoped, reviewed lessons from reproduced results, with a verification step.

Evidence base: `landscape/registry-data-access.md`, `sources/catalog.yaml`, and the first-contact extraction of the ECFSPR 2022 Annual Data Report, all dated 2026-10-04.

---

## 1. `security-browsing` section 4, rule 5: replace an anecdote with tested results

**Why.** The current text says an earlier script got `403` from CFF while an approved fetch tool retrieved some pages. That is true but leaves an agent guessing what to do. We now have a tested table and a machine-readable permission record.

**Insert after the existing rule 5 paragraph:**

> **Tested 2026-10-04 from a datacentre IP with a plain Python client** (evidence and dates in `research/landscape/registry-data-access.md`):
>
> | Site | Result | What it means for agents |
> | --- | --- | --- |
> | `www.cff.org` | `403` on the site root, on `/media/<id>/download` report paths, **and on `/robots.txt`** | We cannot read their crawl policy, so no permission can be inferred. Report PDFs are downloaded by a human in a browser. Never use Playwright or a headless browser to make the request look human; that is the "headless trick" this rule forbids |
> | `databases.lovd.nl` | `402` with a body stating that scraping is not permitted and permission must be requested via their contact form | Explicit refusal. Link only; never fetch |
> | `pr.ecfs.eu` | `200`; `robots.txt` is `User-agent: *` / `Disallow:` (empty) | Automated retrieval of ECFSPR annual report PDFs is permitted. Plain HTTP suffices; browser automation adds nothing |
> | `cftr2.org` | `200`, terms not yet read | Treat as manual until someone reads the terms and records the outcome |
>
> The machine-readable record of these permissions is `research/sources/catalog.yaml`. Its `access` field is a permission, not a hint: an agent may automate retrieval only when it reads `fetch` or `api`. A source absent from the catalog has no permission. See `research/sources/README.md`.

**Verification.** `python research/tools/sources/fetch_sources.py --manual` lists only `access: manual` entries (the CFF reports and CFTR2) and nothing from ECFSPR; `fetch_sources.py` without flags downloads only ECFSPR files.

---

## 2. `cf-research-context`, "Data and community sources": point at the catalog

**Why.** The section lists sources by URL with no access policy, and lists ECFS without a link. Agents reading it will not know that CFTR2's terms are unread or that LOVD refuses automation.

**Replace the bullet list intro sentence with:**

> The workspace's public materials should link to primary sources, not to our own summaries. **Before fetching anything, check `research/sources/catalog.yaml`**: it records, per source, whether automated retrieval is permitted (`fetch`/`api`), must be done by a human (`manual`), needs an application (`request`), or has been refused (`forbidden`). If a source is not in the catalog, add it with evidence before using it.

**Add to the bullet list:**

> - **ECFS Patient Registry (ECFSPR)** — annual data reports back to 2003 as PDFs; <https://pr.ecfs.eu/annual-reports/>. Automated retrieval permitted by their `robots.txt` (checked 2026-10-04).

---

## 3. `cf-research-context`, new subsection after "Disease basics": registry numbers

**Why.** "Numbers need context" is correct but gives no operational guidance for the most common number source we will use. These are the specific traps observed on first contact with a registry report.

**Insert:**

> ### Registry numbers
>
> Registry annual reports are the lab's primary source for population figures. Rules when extracting or quoting one:
>
> - **The denominator changes per indicator, inside one report.** The ECFSPR 2022 report's summary table footnotes three different populations for adjacent rows: everyone registered, only people seen by clinical staff during the year, and only people alive on 31 December. Copy the printed denominator verbatim next to every value. *(Observed 2026-10-04, ECFSPR 2022 Annual Data Report p. 8.)*
> - **Data year comes from the title page, never from the filename or URL.** A CFF report served at a path containing `2019-...` is titled 2020. *(Observed in a search index 2026-10-04; `sources/catalog.yaml` entry `cffpr-adr-2020`.)*
> - **Do not compare across years or registries without a comparability note.** Lung-function reference equations, age bands, inclusion of screen-positive inconclusive diagnoses, and survival methodology can change between reports. The CFF Technical Supplement is the authority for CFF method changes. Cross-registry arithmetic is out of scope until a registry scientist reviews it. *(Inference from registry methods documentation; label `[unverified]` until a specific change is cited.)*
> - **"Median predicted survival" is a model output under stated assumptions, not an observed lifespan.** Quote it with the report's own definition or not at all. *(Standard registry methodology; cite the report's methods section when used.)*
> - **Extracted text is not clean text.** The ECFSPR 2022 PDF's font corrupts `ti`/`tt` ligatures, so labels extract as `Cys-c Fibrosis` and `Introduc?on` while numerals are intact. A local model asked to "tidy" such labels will do so from memory. Normalise deterministically, then have a human confirm the label; never let a model repair a value. *(Observed 2026-10-04.)*
> - **Medians arrive bundled with their interquartile range in one cell**, often split across line breaks. Decide the schema (`value_low`/`value_high` or `value_type: median_iqr`) before extracting. *(Observed 2026-10-04.)*

---

## 4. `ai-loop-council`, "Task packet template": a packet shape for extraction work

**Why.** A lot of the registry work will be delegated to the local worker. The existing packet contract is general; extraction has one specific failure mode worth naming: the model is good at mapping free text to a label and bad at producing a number.

**Add under the template, as a worked variant:**

> **Variant: value extraction from a document.** The worker never produces a number. The deterministic extractor produces candidate `(page, raw_text)` spans; the worker's only job is to map each span to an `indicator_id` from a fixed list, or return `out_of_scope`. The packet contains: the indicator list with definitions, the raw spans inside `<untrusted_page>` tags, and the instruction that abstaining is a correct answer. The verifier checks that every returned id is in the list, that no span was altered, and that abstention rate is reported. A worker that "fixes" a span has failed the task.

**Verification.** The frozen 60-claim set described in the cf-registry-extract scoping memo (a private session document, not in this repo) doubles as the acceptance test: ≥90% correct mappings and zero confident verdicts on the 10 out-of-scope claims.

---

## 5. `model-onboarding`, section 4: a non-clinical CF fixture for the distill benchmark

**Why.** The lab benchmarks local models on non-clinical fixtures. Published aggregate registry indicators are public, aggregate, and have unambiguous right answers, which makes them a better fixture than synthetic text for the "map claim to source, or abstain" task — and it is exactly the task the registry project needs the worker to do well.

**Proposal, design only.** Build `research/tools/sources/fixtures/claim-mapping-v0.json` from the extracted ECFSPR summary indicators: 30 claims the table supports, 20 it contradicts, 10 it does not cover. Score models on verdict accuracy and on abstention discipline separately; a model that answers confidently on the 10 out-of-scope claims fails regardless of its accuracy on the rest. Report with the same small-sample caveats the injection probe uses.

This should not be built until the extraction schema is frozen; otherwise the fixture encodes a schema that will change.

---

## 5b. `security-browsing` section 4: name the API-versus-robots distinction and the enforced rules

**Why.** Section 4 says to check `robots.txt` and read terms, which is right for pages and wrong-by-omission for APIs: NCBI E-utilities' robots file disallows everything and ClinicalTrials.gov's disallows `/api/`, while both publishers document those APIs for programmatic use. An agent following section 4 literally would refuse PubMed. The distinction, and the numeric rules a worker cannot be trusted to infer, are now enforced in code in `skills/cf-evidence-loop`.

**Insert after rule 4.6:**

> 7. **APIs and pages are governed by different documents.** `robots.txt` addresses crawlers of pages. A published API's usage policy addresses API clients, and it is the one that applies to API calls, even where the host's `robots.txt` forbids the path (observed 2026-10-04: `eutils.ncbi.nlm.nih.gov` disallows `/`; `clinicaltrials.gov` disallows `/api/`). Before calling any API, record its terms page and our own rate ceiling in the source catalog; no record, no call. For documents and pages, obey `robots.txt` strictly, and treat a `403` on the robots file itself as a refusal.
> 8. **Numbers, not adjectives.** Default one request per second per source; never more than three per second from this machine to any one host; at most 200 request attempts per process; one honoured `Retry-After` of at most 30 seconds, then stop; any `402`, `403`, or repeated `429` puts that host off limits for the rest of the process; scheduled jobs never more often than hourly per source and never in parallel. These are implemented and tested in `skills/cf-evidence-loop/references/NETWORK-RULES.md`; a tool that cannot demonstrate equivalent rules does not share this IP.

**Verification.** `python -m pytest -q tests` in `skills/cf-evidence-loop` reports the conduct tests passing; `python scripts/evidence_cli.py sources` lists only sources whose catalog entries carry `terms_url` and `max_rps`.

## 6. Housekeeping observed during review

- `AGENTS.md` described `../research` as "multi-agent systems, MCP servers"; updated to match what the repo now holds. *(Applied, since `AGENTS.md` is not under the protected path.)*
- `research/cf-ai-ideas.md` still said the repo was Markdown-only; updated. *(Applied.)*
- `hq/cf-projects/CF-Project-Ideas.md` listed ideas that conflict with `cf-research-context` (dosing tools, individual-level data platforms, a genotype-based drug response predictor). Added a scope-triage table at the top naming each conflict and a salvageable form. *(Applied.)*
- `.github/loop-orchestrator/loop.py` is documented as non-functional and reporting success unconditionally. Nothing here changes that; it is worth either deleting or gating behind a hard failure so it cannot be mistaken for a working verifier.

---

## What was not done, and why

- No skill under `hq/.claude/` was edited. The write gate is correct and this proposal is the intended route.
- No fixture was built. See item 5.
- No model was run. The local GPU is not reachable from the sandbox this was written in; the proposals above are shaped so the lab can run them locally.
