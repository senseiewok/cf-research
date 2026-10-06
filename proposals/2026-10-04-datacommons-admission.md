# Data Commons MCP server: admission review

Status: steps 1 to 4 are done (2026-10-04). Recommendation for the human reviewer: admit the server **for general statistics only** (population and other denominators); it has no cystic fibrosis or sickle cell indicator. The decision is the human's (step 5). Owner repo: research (hq board task T-0048).
Provenance: steps 1 to 4 performed by `claude-sonnet-5-5`; every fact below was read from the cited file, page or saved server response on 2026-10-04. Review tier: full (a third-party server with a key and a network connection).

## Verdicts so far

| Step | Result | Basis |
| --- | --- | --- |
| 1 Provenance | pass | below |
| 2 Terms and rate limits | pass with a gap: **no Data Commons rate limit is published** | below |
| 3 Key handling | pass; the key was supplied only through the git-ignored `.env` and the launcher | below |
| 4 Runtime review (stdio, four queries) | **pass with two disclosed deviations** | below |
| 5 Decision | pending: the human decides | below |

## 1 Provenance

- **Package:** `datacommons-mcp` 1.4.0 on PyPI, uploaded 2026-09-03. Both files (wheel and source archive) were downloaded from the URLs in PyPI's JSON record and their SHA-256 hashes matched the record. Nothing was installed by `pip` until the files had been read.
- **Source:** the declared homepage is `github.com/datacommonsorg/agent-toolkit`, in the `datacommonsorg` organisation. The repository itself was not cloned or read; this review is of the published package.
- **Licence:** Apache-2.0 (full text in the package).
- **Install scripts:** none. The source archive has no `setup.py`; `pyproject.toml` uses the standard `setuptools` backend with no custom hooks; the wheel is pure Python (`py3-none-any`) and cannot run code at install time. The only entry point is `datacommons-mcp = datacommons_mcp.cli:cli`.
- **Runtime dependencies (8 direct):** `fastapi`, `uvicorn`, `fastmcp==3.4.2`, `requests`, `pydantic`, `pydantic-settings`, `python-dotenv`, `google-cloud-storage`. The install into a fresh virtual environment, from binary wheels only so no dependency could run a build script, pulled in 86 packages in total. `pip show datacommons-mcp` reports 1.4.0.
- **Why `google-cloud-storage`:** it is imported lazily, only when the optional setting `DC_INSTRUCTIONS_DIR` starts with `gs://`, to read custom instruction files from a bucket using the machine's default Google credentials. With that setting unset (the lab's choice), it is never used. Keep it unset.
- **What the code does (static read of the 1.4.0 source):** the only hosts in the code are `api.datacommons.org` (API root `https://api.datacommons.org/v2`) and `apikeys.datacommons.org` (a pointer in an error message). A search for process spawning, dynamic code execution, deserialisation, file writes and telemetry or analytics libraries found nothing. Environment variables read: `DC_API_KEY`, `DC_API_KEY_VALIDATION_ROOT`, `DC_AGENT_API_ROOT`, `DC_SEARCH_SCOPE`, `DC_INSTRUCTIONS_DIR`.
- **Two behaviours to plan around.** (a) At startup, `serve` validates the key with one request to the validation root (`api.datacommons.org` by default); that request counts toward our conduct budget. (b) The server calls `load_dotenv(find_dotenv(usecwd=True))` and its settings read `.env`: it loads a `.env` from the **current directory** by itself. Run it from a neutral folder with the key supplied by the launcher, never from a folder that holds a stray `.env`.
- **Instructions the server serves to the agent:** about 4,000 words (a server description, six tool descriptions and three research "skills"). They are third-party instructions and were read before any attachment. They ask the agent to read the skill resources first, never guess place or variable identifiers, and attribute every statistic to its original source. They ask for nothing outside the server's own tools. Untrusted by default all the same.
- **Not verified:** the contents of the GitHub repository beyond the published package, and the behaviour of the hosted endpoint (`api.datacommons.org/mcp`), which the lab does not use.

## 2 Terms and rate limits

Two documentation pages were fetched once each, a minute apart: the API documentation and the MCP server documentation.

- **Keys:** the documentation says an API key is required to authenticate and authorize requests, including MCP server requests at `api.datacommons.org`. Keys are self-service at `apikeys.datacommons.org`, requested for the hostname.
- **Terms:** the footer links Google's general "Terms and Conditions" (`policies.google.com/terms`) and privacy policy. No Data Commons-specific terms page was found on either page.
- **Attribution:** both pages link "Cite Data Commons" at `/attribution.html` on the documentation site. That page was not fetched.
- **Rate limits and quotas:** **not stated on either page.** The catalog records `max_rps: 1.0` as our own ceiling and `published_limit` says plainly that none was published.
- **Licences of the statistics:** Data Commons aggregates many sources, each with its own licence; the per-dataset list is at `datacommons.org/datasets`, which was not read. What the lab needs is each statistic's provenance, which the tool output carries; no blanket licence statement was found.
- **Catalog:** `datacommons` was added to `sources/catalog.yaml` with `access: api`, `auth: api_key_required`, `max_rps: 1.0`, `terms_url`, `claim_label: unverified`, and **no `base_url`**. The evidence skill reaches only sources with a `base_url`, so it cannot reach this one until a human admits it.

## 3 Key handling

- The key lives only in the single lab `.env` (git-ignored; template `.env.example` holds a placeholder). It is loaded into one command by `security-git/scripts/run-with-env.ps1`, which prints nothing about values.
- The agent settings deny reading `.env`; the reviewer never read the file. The key was added by the human and reached the server only as a process variable set by the launcher.
- Acceptance check: a `git grep` for `DC_API_KEY` over the tracked and untracked-but-not-ignored files of the three repos must match only documentation and the placeholder line, never a value.

## 4 Runtime review (run 2026-10-04, 21:58Z to 22:02Z)

The plan was to attach the server to **one** agent client over **stdio only**, with fetch approvals on. Example for Claude Code (adjust the paths to this checkout and the review environment):

```text
claude mcp add datacommons --scope local -- pwsh -NoProfile -Command "Set-Location $env:TEMP; & '<hq>/.claude/skills/security-git/scripts/run-with-env.ps1' -- '<review-venv>/Scripts/datacommons-mcp.exe' serve stdio"
```

Four queries, one minute apart, source and date fields recorded verbatim from the saved responses:

| # | Query (UTC start) | Expected (to be checked, not assumed) | Response |
| --- | --- | --- | --- |
| a | population of the United States, latest (21:58:28Z search, 21:58:32Z observation) | a figure with a US Census source and a date | `Count_Person` for `country/USA`: date `2025`, value `341784857`; `measurementMethod` `CensusPEPSurvey`, `observationPeriod` `P1Y`, `provenanceUrl` `https://www.census.gov/programs-surveys/popest.html`; 14 alternative sources listed (ACS 1-year and 5-year, Wikidata and others) |
| b | population of Washington state, latest (21:59:32Z search, 21:59:35Z observation) | same | `Count_Person` for `geoId/53`: date `2025`, value `8001020`; same primary source and URL; 10 alternative sources |
| c | indicator search: "cystic fibrosis" (22:00:35Z) | nothing CF-specific | none. Topics: Asthma, Tuberculosis, Kidney Disease, Cancer, Diabetes, Hiv, Cyclones, Health, Carbon Monoxide, Dengue. Variables: cyclosporiasis cases, infant respiratory deaths, COPD proportion |
| d | indicator search: "sickle cell disease" (22:01:27Z) | nothing SCD-specific | none. Topics: Kidney Disease, Hiv, Malaria, Tuberculosis, Heart Disease, Cancer, Diabetes, Arthritis, Stroke, Dengue. Variables include deaths from the ICD chapter on diseases of the blood and blood-forming organs and immune disorders (a chapter total, not sickle cell) |

Findings:

- Observations are population-level counts with a named Census method and a date; the source is in the response, so a provider can carry it into an evidence record's `limitations`.
- Neither a CF nor an SCD indicator exists. The nearest hit for (d) is an aggregate blood-disorder death count; it must not be cited as sickle cell data. The registry project's scope does not change.
- Nothing individual-level appeared. The stop pattern matched the word "Patient" inside aggregate variable names in (c); that was read and found to be aggregate counts, so the pattern was narrowed for (d).
- Raw responses are kept outside the repositories (reviewer scratch folder) and are not part of this change.

Deviations, stated plainly:

1. **Scripted client, not an interactive agent client.** The queries were sent by a small throwaway script speaking MCP stdio JSON-RPC to the server through the launcher, single-threaded, no retries, from a neutral working folder. (a) and (b) each used two calls (an indicator search to resolve the DCIDs, then the observation), so the server saw six tool calls in total, plus one key-validation request at each of the two server startups.
2. **Pacing slip on (d).** Query (d) started 52 seconds after (c), not 60: the wait check used a stale threshold. No error, refusal or rate-limit response came back from any call. Recorded here as a conduct slip, not repeated.

No 429, 403 or cooldown response occurred.

## 5 Decision and next step

Step 4 passed. Recommendation: admit for general statistics (population and similar denominators) only, never as a source of CF or SCD figures, and keep the rate limit unpublished-so-be-conservative (1 request per second in the catalog). The decision is the human's. If admitted, the next task is a `datacommons` provider in `cf-evidence-loop` that follows its network rules and returns evidence records whose `limitations` name the underlying source and vintage. That provider is not written here, and nothing in the evidence skill changes until a human says so.

Escalate and stop if: a later install runs a script that was not read; the key cannot be kept out of the repositories; a query returns anything that looks like individual-level data; or the terms restrict redistribution of the statistics in a way that affects the registry table.
