# Small AI tools for CF research

**Proposal inventory, October 4, 2026.** These are hypotheses about useful work,
not validated community needs or clinical products. Check existing alternatives
and seek qualified feedback before choosing a project. Research, not medical advice.

## Three bounded starting points

| Proposal | Smallest useful artifact | Evaluation | Not in scope |
| --- | --- | --- | --- |
| Literature evidence cards | A Markdown schema for study question, population, endpoint, result, source excerpt, date, and limitations | Independently check each field against a small, frozen set of lawful public sources; count unsupported claims and omissions | Treatment recommendations or patient-specific interpretation |
| Historical claim ledger | A sourced timeline that keeps the exact approved statement beside its supporting excerpt | Verify every changed claim and date; ensure educational rewrites preserve meaning | Individual prognosis or unsourced survival statistics |
| Research-tool catalog | A small comparison of existing public CF tooling, data access terms, licenses, and maintenance | Confirm each entry against its source; separate discovery from scientific validation | Collecting patient data or claiming endorsement |

## Proposed first experiment

Start with evidence cards, subject to confirming a real unmet need. Use a small
set of openly accessible sources whose terms permit the intended use. Freeze the
extraction fields and independent reference answers before comparing direct
local work with optional frontier-model planning and bounded local execution.

Measure accepted cards, unsupported statements, missing limitations, total model
tokens, and end-to-end time. A second model's agreement is not a reference answer.
Record negative results; stop or rescope if review effort exceeds the useful work
or if existing tools already meet the need.

## Boundaries and next questions

- Use public aggregate information or synthetic fixtures only. No identifiable health data in repositories, logs, or prompts.
- Notes stay toolchain-free. Code lives under `tools/`, one directory per tool with its own dependency list, README, and tests; see the README conventions. (Superseded the earlier Markdown-only rule on 2026-10-04.)
- Scientific and community reviewers should assess usefulness and interpretation; do not ask contributors to disclose personal medical information.
- Treat [CFF priorities](cff-goals.md) and the [existing-tool landscape](landscape/medical-ai-landscape.md) as starting evidence, not proof that these proposals are needed.

Before building: who would use the artifact, what existing alternative fails,
which sources can lawfully be reused, and what result would make us stop?
