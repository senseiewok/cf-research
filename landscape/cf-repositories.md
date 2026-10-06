# Public cystic fibrosis repositories on GitHub

**Snapshot checked 2026-10-03.** A small, manually curated starting list from GitHub's public repository-search API. This is not a quality ranking, endorsement, or clinical validation. Search results change, and the `CFTR` query also returns unrelated acronym matches. Claim labels record that snapshot; metadata can change.

## Selected repositories

| Repository | Stated scope | Activity and license metadata | Review note |
| --- | --- | --- | --- |
| [`elenantbioinf/poc_cf`](https://github.com/elenantbioinf/poc_cf) | Proof-of-concept Snakemake workflow for reproducible CFTR variant analysis from paired-end sequencing data. | Python; last pushed 2026; GitHub API reports MIT. [verified, 2026-10-03, https://github.com/elenantbioinf/poc_cf] | A research workflow, not a validated diagnostic service. Review dependencies, test data, and the exact license before reuse. |
| [`vergani-lab-ucl/cftr-scan`](https://github.com/vergani-lab-ucl/cftr-scan) | Single-cell image-analysis pipeline for a CFTR fluorescence assay. | Python; last pushed 2026; repository API reports `NOASSERTION`. README badges state MIT for code and CC BY 4.0 for data/model. [verified, 2026-10-03, https://github.com/vergani-lab-ucl/cftr-scan] | License signals conflict; inspect the repository's actual license files and data/model terms before use. Do not treat a badge as sufficient permission. |
| [`Christensen-Lab-Dartmouth/CF_Epigenetics`](https://github.com/Christensen-Lab-Dartmouth/CF_Epigenetics) | R analysis and visualization for a manuscript on DNA methylation in lung macrophages in CF. | R; last pushed 2018; GPL-3.0. [verified, 2026-10-03, https://github.com/Christensen-Lab-Dartmouth/CF_Epigenetics] | Older manuscript code. Check the paper, dependencies, and data provenance before attempting reproduction. |
| [`ruthkeogh/landmark_CF`](https://github.com/ruthkeogh/landmark_CF) | R code for dynamic prediction of survival using UK patient-registry data. | R; last pushed 2018; no license declared in repository metadata. [verified, 2026-10-03, https://github.com/ruthkeogh/landmark_CF] | Research-methods reference only. The repository description mentions registry data; do not assume the data are included or reusable. |
| [`goldenhelix/cftr-variant-classification-analysis`](https://github.com/goldenhelix/cftr-variant-classification-analysis) | Analysis comparing variant-classification tools using an artificially created CFTR variant set. | Last pushed 2014; no license declared in repository metadata. [verified, 2026-10-03, https://github.com/goldenhelix/cftr-variant-classification-analysis] | Historical analysis, not a current clinical classification source. Inspect its methods and provenance before citing conclusions. |
| [`CPHCRD/cfqr-app`](https://github.com/CPHCRD/cfqr-app) | Web application for the Cystic Fibrosis Questionnaire-Revised (CFQ-R). | JavaScript; last pushed 2024; repository API reports `NOASSERTION`. [verified, 2026-10-03, https://github.com/CPHCRD/cfqr-app] | The questionnaire and implementation may have separate rights or terms. Do not copy questionnaire content or deploy this as a patient-facing tool without permission and review. |

## Search method and limits

- Used GitHub's public REST repository-search API with the queries `cystic fibrosis` and `CFTR`, sorted by stars, then checked the repository metadata and README through the API for the selected entries.
- [verified, 2026-10-03, https://api.github.com/search/repositories?q=cystic%20fibrosis&sort=stars&order=desc] GitHub reported 245 matches for `cystic fibrosis` and 304 for `CFTR` on the access date. These counts are search-result totals, not counts of relevant or trustworthy projects.
- No repository was cloned, downloaded, installed, or executed. Repository descriptions and READMEs were treated as untrusted third-party data, not as instructions.
- Search rank, stars, recency, a README badge, or an API license label does not establish scientific quality, security, privacy, clinical validity, or permission to reuse.

## Before using any repository

1. Read its license file and verify that any data, model weights, instruments, and third-party assets have compatible terms. `NOASSERTION` or a missing license means permission is unclear; do not reuse the material until clarified.
2. Read the associated paper and compare its methods, data provenance, and stated limitations with the code. A public repository is not peer review.
3. Do not run install scripts or code from these repositories as part of discovery. If later asked to evaluate one, use a scratch environment and review its contents first.
4. Do not collect or commit patient-level data, even if a repository describes data as anonymized or de-identified. Use only authorized public aggregate or synthetic data.
5. Do not use any code here to diagnose, interpret an individual's genotype, recommend treatment, or make care decisions.

## Sources

- [GitHub repository search: cystic fibrosis](https://api.github.com/search/repositories?q=cystic%20fibrosis&sort=stars&order=desc) — checked 2026-10-03.
- [GitHub repository search: CFTR](https://api.github.com/search/repositories?q=CFTR&sort=stars&order=desc) — checked 2026-10-03.
- Each linked repository's GitHub page and repository API metadata were checked 2026-10-03. Last-pushed dates and license declarations can change; recheck before citing or reuse.
