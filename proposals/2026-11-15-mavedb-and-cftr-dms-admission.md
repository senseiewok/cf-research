# Proposal: Source Catalog Admission for MaveDB (CC0 Zenodo Dump) and CFTR Deep Mutational Scanning Deposits

Status: proposed (written 2026-11-15; awaiting human maintainer decision)

## Summary

This proposal reviews and prepares catalog entries for three major variant-level functional dataset deposits to support the open theratype benchmark (`D2` on the research roadmap):
1. **MaveDB February 2026 Zenodo dump** (CC0 1.0 Universal).
2. **McKee et al. 2023 CFTR deep mutational scanning deposit** (129 variants; *Cell Chemical Biology*).
3. **Tedman et al. 2025 CFTR deep mutational scanning deposit** (232 variants; *eLife* and Mendeley Data).

In Cystic Fibrosis research, patient-level registry data is restricted behind formal data access applications. Open variant-level datasets of deep mutational scanning and pharmacological rescue provide the primary empirical ground needed to benchmark computational prediction models (such as AlphaMissense or language models) against actual CFTR plasma membrane expression and drug responsiveness.

All three sources are proposed under `access: manual` until a maintainer reviews the terms, confirms storage/bandwidth feasibility, and admits them for automated tooling.

---

## 1. MaveDB Bulk Download (Zenodo Dump)

- **Platform & Host:** Zenodo.
- **Concept DOI:** `10.5281/zenodo.11201736` (resolves to the most recent version).
- **Version Issued:** Version 4, issued 2026-02-06.
- **Version DOI:** `10.5281/zenodo.18511521`.
- **Title:** MaveDB Bulk Download.
- **Publisher / Attribution:** MaveDB authors and contributors (hosted on Zenodo).
- **Creators:** MaveDB contributors.
- **Landing Page:** `https://zenodo.org/records/18511521` (and `https://doi.org/10.5281/zenodo.11201736`).
- **Documentation:** `https://www.mavedb.org/docs/mavedb/download.html`.
- **Licence (as stated):** Creative Commons Zero v1.0 Universal (`CC0 1.0`), confirmed via DataCite rights record (`rightsIdentifier: cc0-1.0`).
- **Contents:** Complete database dump containing all published data in MaveDB licensed under CC0:
  - Structured metadata (`main.json`) for experiment sets, experiments, and score sets.
  - Variant score and count files (`.csv`) indexed by score set URNs.
  - Annotation files containing mapped HGVS, VRS identifiers, VEP consequences, gnomAD frequencies, and ClinVar classifications.
- **Relevance to CF Lab:** Contains published multiplexed assays of variant effect (MAVEs), serving as the foundational repository for deep mutational scanning assays.
- **Proposed Catalog Mode:** `access: manual`. The archive is updated quarterly; a human maintainer downloads the bundle or specific CFTR target score sets.

---

## 2. CFTR Deep Mutational Scanning Corrector Response (McKee et al. 2023)

- **Primary Publication:**
  - *Title:* General trends in the effects of VX-661 and VX-445 on the plasma membrane expression of clinical CFTR variants.
  - *Journal:* *Cell Chemical Biology*, 30(6):632–642.e5 (June 2023).
  - *DOI:* `10.1016/j.chembiol.2023.05.001`.
  - *PMID:* 37253358.
  - *PMCID:* PMC10330547.
  - *Authors:* Andrew G. McKee, Eric F. McDonald, W. Dylan Penn, Cortney P. Kuntz, Kevin Noguera, Lauren M. Chamness, Frank J. Roushar, Jens Meiler, Katherine E. Oliver, Lars Plate, Jonathan P. Schlebach.
- **Scope & Dataset Content:**
  - Evaluates plasma membrane expression (PME) across 129 clinical CFTR variants (including 45 previously unclassified variants) using deep mutational scanning (Sort-seq).
  - Measures baseline expression (DMSO vehicle) and rescue by single correctors (VX-661 / tezacaftor, VX-445 / elexacaftor) and combination therapy (VX-661 + VX-445).
- **Deposits & Locations:**
  - *Raw Sequencing:* NCBI BioProject `PRJNA896578` (SRA).
  - *Analysis Code:* GitHub `https://github.com/schlebachlab/RP-dms-analysis` archived at Zenodo DOI `10.5281/zenodo.7799478` under GNU General Public License v3.0 (`GPL-3.0`).
  - *Variant Scores Table:* Published Supplementary Tables (Excel format) in *Cell Chemical Biology* supporting files.
- **Licence (as stated):**
  - Journal article and supplementary tables: Elsevier copyright (non-open access with 12-month NIH Public Access PMC embargo release).
  - Raw sequencing: Public domain (NCBI SRA data terms).
  - Code: GPL-3.0.
- **Proposed Catalog Mode:** `access: manual` (manual download of supplementary data tables).

---

## 3. CFTR Deep Mutational Scanning Calnexin Dependency (Tedman et al. 2025)

- **Primary Publication:**
  - *Title:* General trends in the calnexin-dependent expression and pharmacological rescue of clinical CFTR variants.
  - *Journal:* *eLife*, 14:RP107180 (2025).
  - *DOI:* `10.7554/eLife.107180`.
  - *Authors:* Austin Tedman, John A. Olson III, Minsoo Kim, Catherine C. Foye, et al., Lars Plate, Jonathan P. Schlebach.
- **Scope & Dataset Content:**
  - Deep mutational scanning survey of 232 clinical CFTR variants measuring how ER chaperone calnexin (CANX) regulates plasma membrane expression and pharmacological corrector rescue.
- **Deposits & Locations:**
  - *Processed Quantitative Scores & Source Files:* Mendeley Data, DOI `10.17632/ksfmxv8zts.1`.
  - *Raw Sequencing:* NCBI BioProject `PRJNA1244365` (SRA).
  - *Proteomics:* PRIDE archive `PXD066492` (DIA) and `PXD061664` (DDA).
  - *Analysis Code:* GitHub `https://github.com/schlebachlab/RP-dms-analysis` (`GPL-3.0`).
- **Licence (as stated):**
  - Mendeley Data repository (`10.17632/ksfmxv8zts.1`): Creative Commons Attribution 4.0 International (`CC BY 4.0`), confirmed via DataCite rights record (`rightsIdentifier: cc-by-4.0`).
  - Primary paper (*eLife*): Creative Commons Attribution 4.0 International (`CC BY 4.0`).
- **Proposed Catalog Mode:** `access: manual` (manual browser download from Mendeley Data until a maintainer ratifies).

---

## 4. Scite Retraction and Integrity Verification

All primary DOIs and dataset identifiers were verified live via `evidence_cli.py scite-batch`:
- `10.1016/j.chembiol.2023.05.001`: supporting=9, contrasting=0, mentioning=91; zero retractions, zero editorial notices.
- `10.7554/eLife.107180`: supporting=0, contrasting=0, mentioning=0; zero retractions, zero editorial notices.
- `10.5281/zenodo.11201736`: zero retractions, zero editorial notices.
- `10.17632/ksfmxv8zts.1`: zero retractions, zero editorial notices.
- `10.5281/zenodo.7799478`: zero retractions, zero editorial notices.

---

## 5. Proposed Catalog Entries

```yaml
  # ---------------------------------------------------------------- variant datasets

  - id: mavedb-dump-zenodo
    title: MaveDB Bulk Download (Zenodo archive)
    publisher: MaveDB authors and contributors (hosted on Zenodo)
    kind: dataset
    access: manual
    landing_page: https://zenodo.org/records/18511521
    terms_url: https://zenodo.org/records/18511521
    reference_doi: 10.5281/zenodo.18511521
    claim_label: unverified
    notes: >-
      Bulk download archive of published variant effect data from Multiplexed
      Assays of Variant Effect (MAVEs). Version 4 issued 2026-02-06 (DOI:
      10.5281/zenodo.18511521; concept DOI: 10.5281/zenodo.11201736). Licensed
      under Creative Commons Zero (CC0 1.0 Universal) per DataCite metadata.
      Access set to manual until a human maintainer reviews the bundle size and terms.

  - id: cftr-dms-mckee-2023
    title: CFTR Deep Mutational Scanning Corrector Response (McKee et al. 2023)
    publisher: Cell Chemical Biology (Elsevier BV)
    kind: dataset
    access: manual
    landing_page: https://doi.org/10.1016/j.chembiol.2023.05.001
    terms_url: https://www.elsevier.com/about/policies-and-standards/copyright
    reference_doi: 10.1016/j.chembiol.2023.05.001
    claim_label: unverified
    notes: >-
      Deep mutational scanning plasma membrane expression data for 129 clinical
      CFTR variants under vehicle, VX-661, VX-445, and dual corrector treatment.
      Published in Cell Chemical Biology 30(6):632-642 (2023). Variant scores
      provided in published supplementary tables. Sequencing data at NCBI BioProject
      PRJNA896578; analysis code on GitHub (schlebachlab/RP-dms-analysis; Zenodo
      DOI: 10.5281/zenodo.7799478 under GPL-3.0). Access set to manual.

  - id: cftr-dms-tedman-2025
    title: CFTR Deep Mutational Scanning Calnexin Dependency (Tedman et al. 2025)
    publisher: Mendeley Data
    kind: dataset
    access: manual
    landing_page: https://doi.org/10.17632/ksfmxv8zts.1
    terms_url: https://www.mendeley.com/terms-and-conditions/
    reference_doi: 10.7554/eLife.107180
    claim_label: unverified
    notes: >-
      Deep mutational scanning dataset of 232 clinical CFTR variants measuring
      calnexin chaperone dependence and pharmacological rescue. Dataset hosted
      on Mendeley Data (DOI: 10.17632/ksfmxv8zts.1) under Creative Commons
      Attribution 4.0 International (CC BY 4.0). Accompanying paper published
      in eLife (DOI: 10.7554/eLife.107180; CC BY 4.0). Raw sequencing at NCBI
      BioProject PRJNA1244365; proteomics at PRIDE PXD066492. Access set to manual.
```

---

## 6. Next Steps

1. **2026-11-23 Maintainer Update (`cftr-dms-tedman-2025`):** Maintainer manually downloaded `cftr-dms-tedman-2025` (`DOI: 10.17632/ksfmxv8zts.1`, CC BY 4.0) locally into git-ignored `sources/downloads/tedman-2025/` (`ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx`, `AGM ALT Library Comparison.xlsx`, `cftr_dict3.tsv`) for offline empirical evaluation in `benchmarks/cftr_dms_expression/` (`access` remains `manual`).
2. **2026-11-23 Maintainer Update (`mavedb-dump-zenodo`):** Maintainer manually downloaded `mavedb-dump.2026062418131.tar.gz` (concept DOI `10.5281/zenodo.11201736`, CC0 1.0, `1,907,736,956` bytes, SHA-256 `7101852e929c1458b8e437949c2ec6ac01a4af6ff1e5740f3a656b557be2bb22`, `asOf: 2026-06-24`) and staged `main.json` (`21,120,012` bytes, SHA-256 `a3350ab7538172fd26cea79b5f10fe1e7f48912cdb971b62c152d8296da18fa8`) in git-ignored `sources/downloads/mavedb-2026-06-24/`: `1,254` experiment sets, `1,992` experiments, and `2,777` score sets with `0` CFTR/`P13569`/`ABCC7` or checked CF modifier score sets (the 41 `"fibrosis"` string hits target `MPL`/`P40238`). `access` remains `manual`.
3. `cftr-dms-mckee-2023` (Elsevier supplementary tables; 103 overlapping variants already cross-checked via Tedman et al. 2025's CC BY 4.0 `AGM ALT Library Comparison.xlsx`) remains `access: manual` pending separate maintainer review and staging.
