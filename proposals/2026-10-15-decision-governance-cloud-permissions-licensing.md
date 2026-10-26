# Decision Record: Governance, Cloud-Model Routing, Permission Letters & Skill Licensing (T-0072)

**Status**: ratified (approved by maintainer 2026-10-15)  
**Date**: 2026-10-15  
**Author**: daily-research-sentinel (CRO Governance Round)  
**Task ID**: `T-0072` (cloud-model use on reports with unread terms, permission requests to registries, and skill licensing)  
**Review Tier**: `full`  
**Answer Tier**: `T2`  
**Reference Proposals**: 
- `proposals/2026-10-05-compliance-review.md` (Findings 4, 5, 10)
- `proposals/2026-10-08-registry-permission-letters.md` (Letters 1 & 3)
- `landscape/registry-data-access.md`

---

## Executive Summary

Task `T-0072` addresses three interconnected governance and compliance questions arising from the lab's multi-agent workflows across third-party Cystic Fibrosis research materials and registry reports.

The task acceptance criteria define three human-gated decisions:
1. Whether sources with unread terms route exclusively to local models, and confirming cloud API data-retention settings;
2. Whether to transmit the prepared permission letters to the CF Foundation and the UK CF Trust;
3. Standardizing the licence boundary notice for skills carrying third-party quotations.

This document frames each decision point with factual context, actionable options, risk trade-offs, and concrete recommended rule adjustments. In accordance with maintainer directives and `tasks/HUMAN-TASKS.md`, autonomous agents cannot make legal or policy determinations for the lab; this document prepares the decision package so the maintainer can execute the decisions with clarity.

---

## Decision Point 1: Cloud-Model Routing for Sources with Unread Terms & Account Data Retention

### 1. Context & Observed Facts
- In October 2026, the maintainer downloaded public registry reports (CFFPR, UK CF Registry, Australian ACFDR, Canadian CCFR, European ECFSPR).
- While CFFPR and the UK CF Registry publish explicit permission and copyright notices, the Canadian CCFR and ECFSPR website terms of use remain unread (`sources/catalog.yaml`).
- Cloud models (e.g. Anthropic Claude, Google Gemini) process queries through external API endpoints. Unless accounts have confirmed zero-data-retention and no-training commitments, sending third-party copyrighted text could risk unauthorized ingestion.

### 2. Decision Options for Maintainer
- [ ] **Option 1A: Strict Local-Only Routing for Unread or Restricted Sources (Recommended)**
  - Any source marked `terms: unread` or carrying a restrictive reproduction clause must be processed **exclusively by local models** (`qwen3.8:27b-64k` via Ollama).
  - No raw text excerpts exceeding 8 words from these sources may be dispatched to cloud provider APIs.
  - Cloud models may only receive derived, non-copyrightable facts (e.g., individual verified numbers with source citations) or synthesized paraphrases.
- [ ] **Option 1B: Cloud Processing with Confirmed Zero-Retention Agreement**
  - Verify and log that the lab's commercial API accounts (Anthropic Commercial API, Google Cloud Vertex / AI Studio) enforce contractual Zero Data Retention (ZDR) and do not train on user payloads.
  - Allow bounded cloud model ingestion of unread-terms reports once ZDR is certified in `security-baseline`.

### 3. Proposed Rule Adjustment
Add to `cf-research-context` and `security-browsing`:
```markdown
### Sources with Unread or Restrictive Terms
If a document or registry report carries unread terms or restricts reproduction:
1. Primary extraction must be performed by a local worker model (e.g. `qwen3.8:27b-64k` on loopback `127.0.0.1`).
2. Cloud delegation packets may contain only single isolated numbers, mathematical calculations, and paraphrased summaries below the 8-word overlap limit.
3. Raw page text must never leave the local machine.
```

---

## Decision Point 2: Registry Permission Requests (CFF and UK CF Trust)

### 1. Context & Observed Facts
- **CF Foundation Patient Registry (CFFPR)**: The annual data reports feature a "Figure Permissions" notice stating that use of charts and data requires contacting the Registry team. Automated retrieval to `cff.org` returns HTTP 403.
- **UK Cystic Fibrosis Trust**: The 2024 Annual Data Report copyright notice explicitly states that content may not be used or reproduced in publications without written permission.
- Draft letters have been prepared in `proposals/2026-10-08-registry-permission-letters.md`:
  - Letter 1: CFF Patient Registry (`cffpr-adr-*`)
  - Letter 3: UK CF Trust (`ukcfr-adr-2024`)
- Both actions are tracked as human-only tasks in `tasks/HUMAN-TASKS.md` Section 2.

### 2. Decision Options for Maintainer
- [ ] **Option 2A: Send Letters 1 & 3 from Official Lab Address (Recommended)**
  - Maintainer reviews, fills official contact details, and sends Letter 1 to CFF and Letter 3 to the UK CF Trust via the lab's contact page.
  - Replies are tracked in `landscape/registry-data-access.md`.
  - Until written permission arrives, the lab maintains its strict policy: link and cite only, no table reproduction, max 8-word quotation limit.
- [ ] **Option 2B: Defer Letters and Maintain Pure Citation / Fair Use Mode**
  - Defer sending outbound emails.
  - Strictly limit all registry citations to standard academic fair-use references (Title, Year, Page, Number) with zero table reproduction and zero quoted phrases over 8 words.

### 3. Proposed Rule Adjustment
Update `landscape/registry-data-access.md`:
```markdown
### Registry Permission Status (Updated 2026-10-15)
- CFFPR: Awaiting maintainer transmission of Letter 1. Operating in read/cite mode only; no table redistribution.
- UK CF Trust: Awaiting maintainer transmission of Letter 3. Operating in read/cite mode only; <= 8-word quote limit enforced by `check_source_overlap.py`.
```

---

## Decision Point 3: CC0 Licensing Header for Skills with Quotations

### 1. Context & Observed Facts
- Agent skills in `cf-skills` and `cf-lab` use YAML frontmatter declaring `license: CC0-1.0`.
- CC0 dedicates the author's copyright to the public domain. However, several skills cite definitions, clinical guidelines, or short verbatim phrases from external authorities (e.g. ECFS, FDA, WHO).
- Legally, the lab cannot dedicate third-party copyrighted text to the public domain. In `proposals/2026-10-05-compliance-review.md` (Finding 4), a dual-licensing clarification was introduced for the ECFSPR skill.

### 2. Decision Options for Maintainer
- [ ] **Option 3A: Standardize Dual-Notice Disclaimer Across All Skills (Recommended)**
  - Maintain `license: CC0-1.0` in frontmatter for compatibility with skill marketplaces and agent platforms.
  - Append an explicit copyright boundary notice in every skill body declaring that the CC0 1.0 Universal dedication covers the lab's original compilation, code, and text, whereas cited third-party terms, excerpts, and trademarks remain under their respective ownership and fair use boundaries.
- [ ] **Option 3B: Complete Elimination of Third-Party Quotations**
  - Enforce zero literal quotations across all skill files: rewrite all definitions into 100% original lab paraphrases.
  - Maintain pure CC0-1.0 without disclaimer blocks.

### 3. Proposed Rule Adjustment
Add to `ai-provider-compatible-skills` and `cf-lab/AGENTS.md`:
```markdown
### Skill Licensing & Third-Party Text
When a skill carries short third-party quotations, include the standard copyright reservation note distinguishing the author's CC0 dedication from quoted third-party text.
```

---

## Maintainer Action Sign-Off

Ratified by maintainer directive on 2026-10-15:

- **Decision 1 (Cloud Routing)**: Approved **Option 1A with Fallback to Option 1B**  
  *Policy*: Local-first extraction via `qwen3.8:27b-64k` loopback. Bounded cloud fallback to Frontier models (Gemini & Claude team) under certified Zero Data Retention (ZDR) and non-training terms.
- **Decision 2 (Permission Letters)**: Approved **Option 2B (Retaining Manual Tasks)**  
  *Policy*: Operate under strict academic citation and <= 8-word quote limits without outbound letters today; retain permission letters in `tasks/HUMAN-TASKS.md` Section 2 as pending human actions for future dispatch.
- **Decision 3 (Skill Licensing)**: Approved **Option 3A (Dual-Notice CC0 Dedication & Reservation)**  
  *Policy*: Standardize dual-notice disclaimer across skills (`license: CC0-1.0` covers original lab code and text; third-party quotes remain under fair use with tight quotation hygiene).

**Maintainer Decision Status**: Ratified  
**Date**: 2026-10-15  
