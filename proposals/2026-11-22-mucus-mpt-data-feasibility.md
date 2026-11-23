# Feasibility & Data-Gate Record: Public Airway-Mucus MPT Data and CFTR DMS Admission Gates (2026-11-22)

Status: review record (written 2026-11-22; Council Item #3 raw MPT benchmark verdict: BLOCKED — could not test, no admitted raw airway-mucus MPT dataset in `sources/catalog.yaml`)

Independent open-source research record, not medical advice. Zero external datasets downloaded and zero catalog permissions modified.

---

## 1. Deterministic Audit of `sources/catalog.yaml` (`355` entries, audited 2026-11-22)

Before running any new mucus-transport benchmark (Council Item #3), we audited all `355` entries in `sources/catalog.yaml` and checked local data artifacts in this repository:

- **Zero admitted raw airway-mucus or CF sputum MPT datasets (`0` eligible):** Only two catalog entries combine `kind: dataset` with `access: fetch` (`natural-earth-110m-countries` geographic boundaries and `clinvar-bulk` variant classifications), and zero `access: api` or `access: fetch` dataset entries provide raw ex vivo human airway-mucus or CF sputum multiple-particle-tracking (`MPT`) trajectories.
- **Mucus literature entries are manual article records only:** Existing mucus-delivery catalog entries (`conte-2022-acsami`, DOI `10.1021/acsami.1c14975`; `comegna-2021-scirep`; and `alton-2024-erj-rsiv-fhn-lentivirus-preclinical`) are `kind: journal_article` with `access: manual`, not open raw particle-tracking trajectory tables.
- **Existing `zenodo5533877_mucus_gel.npz` is out of domain:** The local artifact `tools/curvature-screen/benchmark/out/zenodo5533877_mucus_gel.npz` (Zenodo record `5533877`, Larobina et al. 2021, CC BY 4.0) contains Rheo-DLS intensity autocorrelation (`g₂(τ)−1`) measurements on purified pig gastric mucus gel. As documented in `tools/curvature-screen/BRIEFING.md`, it is purified pig gastric mucus rather than human airway mucus or CF sputum, and it does not record single-particle MPT trajectories.
- **External Zenodo SEM workflow (`16272454`):** Zenodo record `16272454` provides a Fiji/Weka SEM mucus-network imaging workflow with example documentation; it is neither admitted in `sources/catalog.yaml` nor a matched transport-trajectory cohort.

---

## 2. Decision Gate (`BLOCKED / could not test`) & Admission Requirements for Candidate MPT Data

Under the Discovery Protocol (`cf-lab/docs/lab-focus-2026-11-14.md`, Step 4), when `sources/catalog.yaml` admits no compatible measured dataset, the outcome is recorded as **`could not test: no data`** rather than generating synthetic particle tracks or relabeling pig gastric hydrogel Rheo-DLS data as CF airway validation.

Any future candidate public airway-mucus MPT dataset must satisfy all of the following gates before benchmarking:

1. **Catalog & license admission:** Explicit reuse license, terms review, and human maintainer ratification in `sources/catalog.yaml` (`access: fetch` or `api`, `claim_label: verified`, or explicitly approved local download provenance with SHA-256 checksum).
2. **Specimen-level holdout identifiers:** Explicit donor, specimen, or biological-batch IDs so train and held-out evaluations partition whole specimens rather than fragments of the same trajectory or slide.
3. **Physical & imaging metadata:** Particle hydrodynamic diameter, surface chemistry/PEG grafting and measured zeta potential, lag-time (`τ`) schedule, temperature, and localization-error, drift, and immobile-particle censoring metadata.
4. **Tissue & disease scope labeling:** Strict separation between non-CF human airway mucus, cell-free hydrogel simulants, animal mucus, and ex vivo human CF sputum (`non-pooling invariant`).

**Bounded next step for Council Item #3:** Prepare a targeted source-admission review of Kaler, Joyner & Duncan (2022, *APL Bioengineering*) or an equivalent public ex vivo human airway-mucus MPT data deposit to verify raw trajectory availability, metadata completeness, and reuse terms.

---

## 3. Companion Council Item #1 Data Gate (`proposals/2026-11-15-mavedb-and-cftr-dms-admission.md`)

Council Item #1 (quantitative CFTR plasma-membrane expression and pharmacological rescue benchmarking) remains gated by `proposals/2026-11-15-mavedb-and-cftr-dms-admission.md` (`Status: proposed`):

- `mavedb-dump-zenodo` (Zenodo `18511521`, concept `11201736`, stated CC0 1.0): `access: manual`, `claim_label: unverified`.
- `cftr-dms-mckee-2023` (DOI `10.1016/j.chembiol.2023.05.001`, Elsevier supplementary tables): `access: manual`, `claim_label: unverified`.
- `cftr-dms-tedman-2025` (DOI `10.7554/elife.107180`, Mendeley Data `ksfmxv8zts.1`, stated CC BY 4.0): `access: manual`, `claim_label: unverified`.

**Guardrail note:** `benchmarks/theratype/cftr_theratype_benchmark.py` Control #5 (`unadmitted_manual_source_rejected`) intentionally asserts fail-closed `PermissionError` against all three IDs until a human maintainer ratifies them. If a maintainer later ratifies `cftr-dms-tedman-2025` or `mavedb-dump-zenodo`, update Control #5 to test a dedicated synthetic unadmitted catalog entry rather than relaxing the allowlist rule (`access in {'fetch', 'api'}` and `claim_label == 'verified'`).
