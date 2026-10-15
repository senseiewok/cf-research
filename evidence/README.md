# Evidence Packages (`evidence/`)

This directory stores durable **Task Evidence Packages** and **Sentinel Verification Run Manifests**.

## Purpose

When research is performed at task conception, compressing findings into a single-line string in `tasks/board.json` risks losing table numbers, verbatim citations, cohort exclusion filters ($N$), and regulatory dates. When an implementation agent picks up the task days later, it is forced to re-excavate source PDFs, risking token burn, denominator dilution, or hallucinated claims.

An **Evidence Package** anchors primary research artifacts directly to a board task via `tasks/board.json`'s optional `evidence_package` field. Downstream agents can directly read these packets, inject exact numbers into their prompt context, and verify claims deterministically.

## When to Create an Evidence Package

- **Recommended for**: Curative biophysical tasks (Track A), clinical translation & patient platform tasks (Track B), registry epidemiological metrics, regulatory drug label extractions, and multi-source evidence syntheses.
- **Not required for**: Routine code hygiene, dependency bumps, CI configuration, unit test updates, or internal governance.

## File Naming Convention

```text
evidence/<YYYY-MM-DD>-<task_id>-<descriptive-slug>.md
```
Example: `evidence/2026-10-15-t0146-cffpr-percentages.md`

Broad multi-task landscape dossiers live beside this directory in `landscape/` (e.g., `landscape/multi-registry-eligibility-2026.md`).

## Packet Structure & The 7 Mandatory Sections

Every task evidence packet must carry YAML frontmatter and the following 7 standard sections:

```markdown
---
task_id: T-0146
title: US CFFPR 2024 Modulator Ineligibility Demographics
date: 2026-10-15
answer_tier: T2
review_tier: full
catalog_sources:
  - cffpr-adr-2024
claims_file: docs/cffpr-waiting-claims.json
status: verified
---

# Task Evidence Packet: <task_id> - <title>

## 1. Curative Mission & Biophysical Relevance
State whether this work is Track A (Curative Biophysics) or Track B (Clinical Translation). Describe how this evidence informs the search for curative therapies (e.g. PTC readthrough, alternative conductance) or prevents misrepresentation in patient translation.

## 2. Primary Sources & Exact Locators
List every primary source ID from `sources/catalog.yaml`, document date/version, and exact structural locator (table name, figure number, page, section).

## 3. Verbatim Excerpts & Quantitative Extractions
Provide verbatim quotes from the primary source text and a clear table of extracted numerical metrics.

## 4. Denominators, Exclusions, & Cohort Boundaries ($N$)
Specify the exact denominator $N$ for every percentage or sample statistic. Detail who was included, who was excluded, and age/encounter/insurance filters.

## 5. Invariants Against Misuse (The Guardrails)
Explicitly document negative rules:
- Prohibit false population equivalence (e.g., sample percentage != total registry percentage).
- Enforce the Non-Pooling Imperative (never sum or blend stats across different national registries).
- Record temporal boundaries (e.g. drug approvals not active during the reporting calendar year).

## 6. Limitations & Explicit "Not Stated" Register
List what the primary source does **not state** (e.g. unstated genotype sub-distributions, unstated socioeconomic details).

## 7. Deterministic Verification Command
A single command that verifies the claims against primary source evidence (e.g., `check_claims.py` command or test suite).
```

## Anti-Hallucination Invariant

In accordance with `cf-research-context`, every factual claim extracted in an evidence package must trace to an official source catalog entry. Never infer or pool numbers across registries.
