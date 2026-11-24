# Proposal: Single-Track External Replication & Collaboration Brief (`Council Item #7`)

- **Status:** proposed (awaiting maintainer approval; zero outreach sent)
- **Written:** 2026-11-23
- **Release status:** **internally checked, not independently reproduced**
- **Target track:** Single-Track Primary Pilot (`benchmarks/theratype/`, `T-0925` / `PRED-0010`)
- **Privacy boundary:** Zero patient records, zero PHI, zero clinical advice

## 1. Purpose & Why Start With One Track

An outside analyst running a small, self-contained check without private assistance provides stronger validation than additional internal model runs. Rather than asking external collaborators to test three unfinished tracks at once (where manual deep mutational scanning or human sputum particle-tracking files are still awaiting maintainer admission), this brief locks **one primary pilot track** (`benchmarks/theratype/`) that runs completely offline from checked-in open fixtures while transparently marking what is not yet admitted.

## 2. Reproducibility Readiness Inventory

| Asset | Code Available | Frozen Fixture / Audit | Empirical Holdout Admitted | Executable Offline | Role in Outreach |
| --- | --- | --- | --- | --- | --- |
| `benchmarks/theratype/` | Yes | Yes (`2026-11-21.v2`, `exposed_development_reference`) | No (`experimental_rescue_status = unknown`; Item #1 DMS pending) | Yes | **Primary collaboration pilot** |
| `benchmarks/mucus_transport/` | Yes | Yes (`audit_mucus_barrier.py`, 5 planted defects, synthetic identifiability) | No (Item #3 CF sputum MPT pending admission) | Yes | Secondary mathematical/correctness audit |
| `tools/readiness/` | Yes | Yes (`cf_ml_benchmark.json`, `stress_eval_synthetic.json`) | No (dual-reviewer unseen paper set pending) | Yes | Secondary literature-extraction harness |
| Item #1 CFTR DMS (`mavedb-dump-zenodo`, McKee 2023, Tedman 2025) | No (`PermissionError` gate active) | Metadata only | No (`access: manual`) | No | Blocked until maintainer ratification |
| Item #3 CF sputum MPT dataset | Screener only | No admitted human MPT file | No (`BLOCKED — could not test`) | No | Blocked until maintainer dataset decision |

## 3. Independence & Scope Boundary Matrix

Running the locked scripts below verifies software reproducibility, split accounting, and label provenance across independent people; it does **not** constitute biological or wet-lab reproduction.

| Independence Dimension | Status in This Pilot | Meaning for External Reporting |
| --- | --- | --- |
| **People** | `independent` (invited CFTR specialist + outside quantitative analyst) | Different readers outside the lab evaluate endpoints, provenance, and code behavior |
| **Code** | `shared-locked` | External analyst executes our public repository scripts and independent audit script as committed |
| **Data** | `shared-locked` | Open FDA SPL label concordance panel (`30` variants) is shared and already marked `exposed_development_reference` |
| **Laboratory / Assay** | `not-applicable-no-wet-lab` | Zero new cell-line, organoid, or human measurements are generated in this step |

## 4. Primary Pilot Question & Reproduction Protocol (`benchmarks/theratype/`)

### Focused Research Question
Does separating **US FDA Structured Product Label listing concordance** (`listed_on_fda_spl` vs `not_listed_on_fda_spl` across `TRIKAFTA` v23 and `ALYFTREK` v9) from **measured experimental CFTR rescue** (`unknown`), paired with domain-disjoint (`NBD1`, `TMD2`, `NBD2_C_terminal`) and 4-fold Leave-One-Family-Out splits, form a scientifically valid negative baseline showing that sequence pathogenicity scores (`AlphaMissense`) do not substitute for modulator-response labels?

### Baseline & Current Internal Output
- **Pre-registered finding (`PRED-0010`):** Refuted (`AlphaMissense` abstention-aware ROC-AUC `0.5101` on the held-out `TRIKAFTA` v23 test split vs deterministic sequence-consequence baseline `0.6818`; `0.0` experimental rescue coverage).
- **Important limitations:** This fixture is already `exposed_development_reference` and tracks **regulatory label listing**, not cellular `%WT` chloride rescue. Five benchmark variants differ between `TRIKAFTA` v23 and `ALYFTREK` v9 labels, illustrating why `"not_listed_on_fda_spl"` cannot be treated as measured nonresponse.

### Offline Reproduction Commands (Clean Repo Root)

```bash
git clone https://github.com/senseiewok/cf-research.git
cd cf-research
python -m unittest discover -s benchmarks/theratype -v
python benchmarks/theratype/cftr_theratype_benchmark.py --check --falsify
python benchmarks/theratype/independent_audit.py --check
```

## 5. Specific Ask to External Collaborators

1. **To a CFTR biology / theratyping specialist:**
   - Critique the domain boundaries, sequence-consequence baseline, and `label_status` vs `experimental_rescue_status` separation in `benchmarks/theratype/`.
   - Identify any variant-classification error, splicing-boundary oversight, or flaw in the pending DMS admission plan (`proposals/2026-11-15-mavedb-and-cftr-dms-admission.md`) before human data admission unlocks held-out testing.
2. **To an independent quantitative analyst:**
   - Run the commands above cleanly from a fresh checkout with no private assistance and report any environment assumptions, edge-case bugs, or discrepancy in `independent_audit.py` or `benchmarks/mucus_transport/audit_mucus_barrier.py`.
3. **If no collaborator run is completed:** Keep every release explicitly labeled **"internally checked, not independently reproduced."**

## 6. Governance, Community Review & Outreach Gate

- **No outreach sent:** Sending invitations, emails, or external review requests requires explicit maintainer approval first.
- **Community review & compensation:** Include people affected by CF or community reviewers when designing plain-language explanations of benchmark results, and compensate external/community reviewers fairly where feasible rather than requesting unpaid labor.
- **No clean-room distraction:** Consolidating unrelated repository modules is **not** a prerequisite for locking or sharing this single-track reproduction packet.
