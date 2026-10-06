# Survey: MCP servers and APIs worth adding to the lab

Status: proposal, 2026-10-04. Discovery is not endorsement, security review, or validation (`security-browsing` section 9). Every third-party server below must pass `security-runtime` and `model-onboarding` before an agent in this workspace gets it as a tool. Access labels are from reading documentation, not from a lawyer.

## 1. The honest comparison first: BioMCP

BioMCP (GenomOncology, MIT licence) is an MCP server and CLI over biomedical sources. Its documentation describes a single command grammar reaching roughly thirty sources including PubMed/PubTator3, Europe PMC, ClinicalTrials.gov, ClinVar, MyVariant.info, bioRxiv/medRxiv, NIH RePORTER, Drugs@FDA and openFDA labels, with trial, variant, drug and article entities. `[verified: project README and registry listings read 2026-10-04; not installed or run]`

That overlaps most of `cf-evidence-loop`. What our skill does that BioMCP does not claim to:

| Property | `cf-evidence-loop` | BioMCP (as documented) |
| --- | --- | --- |
| Permission gate | Catalog `access` field is the only permission; refusal is a unit test | Reaches its configured sources; no per-source permission record in the lab's terms |
| Record shape | Every answer is an `Evidence` record with a mandatory non-empty `limitations` field | Returns structured results with source attribution |
| Model in the loop | None, by design | Includes a "think" tool for multi-step reasoning |
| Evidence ledger | `--ledger` appends JSONL | Not a stated feature |
| Scope | Eight providers, CF/SCD examples, tiny | Broad, oncology-strong |

Recommendation: **evaluate BioMCP as a discovery tool and keep `cf-evidence-loop` as the evidence-record layer.** Concretely, an agent may use BioMCP to find candidates and must use `cf-evidence-loop` (or a direct API call with the same record shape) to produce the citation that ships. Before any agent gets BioMCP as a tool: run `model-onboarding` section 1 (provenance: pin a release, read the install scripts), the injection probe, and the `security-runtime` MCP checklist; confirm which of its upstream sources need keys and which the lab has catalogued; and read its handling of the OpenAI key variable its install notes mention, since the lab's loop must stay local-first.

## 2. Other MCP servers to consider

| Server | What it adds | Lab fit | Before adopting |
| --- | --- | --- | --- |
| GitHub's official MCP server | Issues, PRs, code search through the API rather than scraping | High: the board, PRs, and the three repos | Scope the token to the lab's repos only; read-only where possible |
| Playwright MCP (Microsoft) | Browser automation as a tool | Medium: UI testing for `.ai` (T-0016) | **Only for the lab's own site.** The catalog gate does not apply to a browser tool; `security-browsing` rule 4.5 does. Never point it at a publisher that refused automated access |
| Filesystem reference server | Scoped file access | Medium | Pin the root to a repo, never a home directory |
| Any PubMed / ClinicalTrials single-purpose server | Narrower than BioMCP | Low: already covered twice | Skip unless BioMCP fails review |

No MCP server is attached to this Claude Science session beyond Gmail, Calendar and Drive, so none of the above was exercised here. `[verified: connector list checked 2026-10-04]`

## 3. APIs not yet in the catalog, by value

| API | What it gives | Access | Why it matters |
| --- | --- | --- | --- |
| **MyVariant.info / MyGene.info** | Aggregated variant and gene annotation (ClinVar, gnomAD frequencies, dbSNP) in one call | public, keyless | One record per variant instead of three providers; aggregate only |
| **gnomAD GraphQL** | Population allele frequencies for CFTR and HBB variants | public; terms to read | Carrier-frequency claims need a population denominator; this is it |
| **Semantic Scholar Graph API** | Citation graph, influential-citation flags, recommendations | public; key recommended | The "paper trail" layer: who built on a result |
| **Unpaywall** | Legal open-access location for any DOI | public; requires a contact email parameter | Turns "cited" into "readable" without touching a paywall |
| **OpenAlex** | Works, concepts, institutions, funders | catalogued; **key required** | Already in the catalog; the provider is unwritten because no key is configured |
| **PubTator3** | Entity annotations (gene, disease, chemical, variant) on PubMed abstracts | public, NCBI | Makes `pubmed` results filterable by what they are about, not what they say |
| **DailyMed SPL API** | The canonical label document, by set_id | public | openFDA already returns the set_id; this fetches the authoritative copy |
| **ISRCTN registry API** | UK and international trials not on ClinicalTrials.gov | public | Non-US trials layer |
| **WHO ICTRP** | Cross-registry trial search | bulk export only; no live API | Layer 3 completeness; batch job, not a provider |
| **EMA medicines data** | EU approvals and EPARs | downloadable tables; terms to read | "Approved in Europe" sentences |
| **Orphanet / Orphadata** | Rare-disease nomenclature, prevalence estimates with sources | public; CC BY | Shared vocabulary across CF, SCD and other rare conditions |
| **NCBI MedGen / MeSH APIs** | Concept identifiers | public, NCBI | Makes indicator slugs interoperable |
| **CFF drug development pipeline** | What is in trials for CF, by phase | **manual**, www.cff.org is 403 | Human-maintained snapshot only |

Not recommended as data sources, for the reasons in `security-browsing`: patient forums and social media, any individual-level dataset, press releases as evidence.

## 4. For the sickle cell community specifically

The evidence loop is disease-agnostic; the SCD-specific sources that would strengthen it are:

| Source | What it gives | Access |
| --- | --- | --- |
| ClinVar and gnomAD for HBB | Variant classification and allele frequency | as above |
| Cochrane Cystic Fibrosis and Genetic Disorders Group | Historically the Cochrane group covering **both** CF and haemoglobinopathies; its SCD reviews reach PubMed and the `reviews` command | via PubMed abstracts today |
| ASH clinical practice guidelines on SCD | The field's recommendations, layer 6 | publisher terms to read |
| NHLBI sickle cell resources; Sickle Cell Disease Association of America | Community and funding context | public pages; `manual` until terms are read |
| Farooq et al., JAMA Netw Open 2020, doi:10.1001/jamanetworkopen.2020.1737 | The funding comparison between SCD and CF | resolved through the tool 2026-10-04: open access, 140 Europe PMC citations, no retraction notice |

## 4b. Paid MCP servers worth the money, and which are not

Prices below were read on 2026-10-04 from vendor documentation where it exists and from third-party pricing indexes otherwise; the indexes disagree with each other on Consensus and Elicit, so treat those two as `[unverified]` until the vendor page is read at purchase time. All three are **hosted** MCP servers: the vendor's servers make the upstream requests, authenticated by OAuth or key, metered by credits. That has two consequences the lab cares about. The operator's IP is not the one making bulk requests, which is the concern NETWORK-RULES.md exists for. And every query leaves the machine to a third party, so `security-browsing` section 2 applies in full: public research questions only, never private material.

| Service | What it adds that nothing free does | Pricing read | Fit |
| --- | --- | --- | --- |
| **Scite** (hosted MCP; OAuth or API key) | **Smart Citations**: each citing statement classified as supporting, contrasting, or mentioning, across a very large citation index. This is the only paid capability that directly serves claim verification: it tells you whether later literature *supported or contradicted* a finding, not just how often it was cited. Also a reference-check API and retraction flags on collections. | Vendor docs: Basic $20/mo (250 MCP credits), Pro $50/mo (2,500 MCP credits, self-service API key), Team $50/seat, Enterprise custom. The Search API for research use needs a separate licence. `[verified: docs.scite.ai pricing page, 2026-10-04]` | **Highest.** One Pro seat is the first paid thing the lab should try. Wire it as a provider: `(doi) -> supporting / contrasting / mentioning counts + the statements`, with the limitation that classification is the vendor's model output |
| **Consensus** (hosted MCP on paid plans) | Yes/no "Consensus Meter" over peer-reviewed papers; citation grounding that shows the quote behind each cited claim. | Indexes report Free; Pro around $10-20/mo with ~250 API+MCP uses; Deep around $45-65/mo with ~1,000 uses. `[unverified: third-party indexes disagree]` | Medium. A synthesis layer: by the lab's own rule a model's summary is not evidence, so it is a discovery tool like BioMCP. The grounding quotes are the useful part |
| **Elicit** (hosted MCP and REST API on Pro) | Systematic-review workflow (screening thousands of papers, structured extraction columns), clinical-trial search, SOC 2. | Indexes report Free tier; Pro around $49/mo billed annually with API and MCP. `[unverified]` | Medium for one specific job: if the lab ever runs a formal review (e.g. the registry-comparability literature), this is the tool. Otherwise it duplicates PubMed plus synthesis |
| Web of Science, Scopus, Dimensions | Curated citation indexes with institutional pricing | Institutional; not for an unaffiliated lab | Low. Europe PMC and Crossref cover the lab's citation needs; OpenAlex covers the rest once a key is configured |

Before paying for any of them: run the same admission review as BioMCP (`security-runtime` MCP checklist, what leaves the machine, what the vendor retains), and write the provider so its records carry a `limitations` string saying that the classification or synthesis is the vendor's model output and must be checked against the cited statement. The lab buys data about citations; it does not buy conclusions.

## 4c. Google, since the lab already pays for Gemini

| Offering | What it is | Cost | Fit for this lab |
| --- | --- | --- | --- |
| **Data Commons MCP Server** (official, open) | Google's open knowledge graph of public statistics (census, health, economics) exposed as an MCP server; natural-language indicator search then observation retrieval, with sources attached. Starter packages on PyPI, quickstarts for Gemini CLI and ADK. `[verified: Google announcement coverage and GitHub, 2026-10-04; not installed]` | Free | **High for layer 5 denominators.** Prevalence and carrier-frequency sentences need a population base; this gives US and state populations with provenance. Does not hold CF or SCD registry data as far as the documentation shows `[unverified]`; it complements the registry extraction, it does not replace it. Admission review still applies, but it is Google-published and open |
| **Gemini Deep Research agent** (Interactions API) | A managed long-horizon research agent that plans, searches, and synthesises a cited report in the background; two versions (standard and Max); a *collaborative planning* mode where the developer reviews and refines the plan before execution; MCP tool support per launch coverage. `[verified: ai.google.dev docs, 2026-10-04]` | Paid Gemini API tiers; usage-based | **High, for exactly one role.** The website already names the lab's research direction: plan with a frontier model, execute locally. Deep Research in collaborative-planning mode *is* the frontier-planning leg, and it browses from Google's infrastructure, not from the lab's IP. By the lab's rules its report is not evidence: every DOI, trial, and approval it cites gets run through `cf-evidence-loop` before anything ships. Discovery at scale, verification at home |
| Gemini API grounding with Google Search | A model tool that cites web results | Included in API pricing | Medium. Same rule: discovery, then verify |
| NotebookLM | Grounded reading and Q&A over documents you upload, in the Gemini subscription | Included | Medium, as a **human** tool: load the CFF Technical Supplement and the ECFSPR methods chapter and ask comparability questions while writing `comparability_notes`. Not for agents; no public API confirmed `[unverified]` |
| **Google Scholar** | Search over scholarly literature | Free to browse | **Never automate.** There is no API, and the terms of service prohibit automated access. A scraper here is the fastest known way to get an IP blocked. Use PubMed, Europe PMC, Crossref, OpenAlex |
| Google Patents public datasets (BigQuery) | Patent full text and metadata | BigQuery pricing | Low until a patent question exists |

The combination worth piloting: Deep Research with collaborative planning for the broad question ("what has changed in CF registry survival methodology since 2019?"), the plan reviewed by a human, the cited report handed to the local loop, and each citation turned into an evidence record. That is the frontier-plans-local-execution experiment the website promises, with a verification step the website does not yet mention.

## 4d. Microsoft, and the other large vendors

Microsoft publishes an official catalog of its MCP servers at `github.com/microsoft/mcp`; the lab develops with GitHub Copilot and VS Code, so these are relevant as **development tooling**, not as research sources. `[verified: catalog and Microsoft developer blog read 2026-10-04]`

| Server | What it is | Fit |
| --- | --- | --- |
| Microsoft Learn Docs MCP (remote, `https://learn.microsoft.com/api/mcp`) | Cloud-hosted semantic search over official Microsoft documentation, returning content chunks with titles and URLs | Medium. Keeps Copilot current on .NET, PowerShell, VS Code and Azure APIs the lab's scripts touch. Hosted by Microsoft; nothing leaves but the query |
| GitHub MCP server (GitHub is Microsoft) | Issues, PRs, code search through the API | High for the board and PR workflow this packet creates. Scope the token to the lab's repositories |
| Playwright MCP (Microsoft) | Browser automation as an agent tool | Medium, **for the lab's own site only** (T-0016 rose count). Never pointed at a publisher; the catalog gate does not govern a browser tool, `security-browsing` rule 4.5 does |
| Azure MCP, Microsoft Fabric MCP, Foundry MCP (remote `https://mcp.ai.azure.com`), Azure DevOps remote MCP | Cloud and infrastructure management | Low. The lab has no Azure estate |

Microsoft has no research-data offering comparable to Data Commons. Microsoft Academic Graph was retired in 2021; its successor is OpenAlex, which is already in the catalog and awaits a key.

**Other vendors, briefly.** AWS publishes official MCP servers for its own services and hosts the Registry of Open Data, which includes public health datasets by download rather than by MCP; relevant only if a specific dataset is needed `[unverified: not read this session]`. Hosted search APIs built for agents (Perplexity, Exa, Tavily and similar) move web crawling to the vendor's infrastructure, which protects the lab's IP, but they are discovery tools under the lab's rules and their results must still be verified through the loop. Nothing from these vendors changes the recommended order below.

## 5. Recommended order

1. BioMCP admission review (T-0042). Decide complement-or-not with evidence, not on reputation.
2. MyVariant.info and gnomAD providers for the variant layer; both are aggregate.
3. Unpaywall and Semantic Scholar for the paper-trail layer; both need a contact or key, which the catalog's `auth` field already models.
4. OpenAlex provider once a key is configured.
5. Everything else as a claim needs it.
