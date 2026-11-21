# Open CFTR Theratype Benchmark (`T-0925` / `PRED-0010`)

**Research benchmark only; not clinical or diagnostic guidance. Zero patient or PHI data.**

Evaluates whether open variant-effect predictors — led by **AlphaMissense** (Cheng et al., *Science* 2023, DOI `10.1126/science.adg7492`; CFTR benchmark by McDonald et al., *PLOS ONE* 2024, DOI `10.1371/journal.pone.0297560`, `PMC10810519`), alongside **REVEL** and **CADD Phred** — predict Cystic Fibrosis Transmembrane Conductance Regulator (`CFTR`) modulator responsiveness across missense, in-frame indel, splicing, and Class I premature-stop/frameshift alleles when compared against a deterministic sequence-consequence class baseline.

## Primary Open Data & Catalog Provenance (`2026-11-21.v2`)

Every value in `theratype_open_annotations_2026-11-21.json` is recorded from admitted `sources/catalog.yaml` APIs (`myvariant`, `openfda-label`, `clinvar`, `ncbi-eutils`, `uniprot-p13569`, `europe-pmc`):

- **Modulator responsiveness ground truth (`openfda-label`)**: US FDA Structured Product Labels for `TRIKAFTA` (NDA 212273/217660, v23, effective `2026-06-22`, Table 6: 521 variants), `ALYFTREK` (NDA 218730, v9, effective `2026-08-26`, Table 5: 564 variants), `KALYDECO` (NDA 203188, v36, effective `2026-03-23`, Table 3: 97 variants), and `SYMDEKO` (NDA 210491, v22, effective `2026-08-21`, Table 6: 154 variants).
- **Variant pathogenicity & conservation scores (`myvariant`)**: Canonical `ENST00000003084` / `P13569` `dbnsfp.alphamissense`, `dbnsfp.revel`, and `cadd.phred`. All non-missense alleles (`F508del`, `I507del`, intronic/canonical splice variants, and Class I nonsense/frameshift alleles) strictly carry `None` (`"not stated"`) for missense-only predictors (`alphamissense_score` and `revel_score`), and uncataloged foundation scores (`evo2_zero_shot_score`, `alphagenome_avi_score`) remain `None`.
- **Fail-closed manual-access gate (`T-0928`)**: `mavedb-dump-zenodo`, `cftr-dms-mckee-2023` (DOI `10.1016/j.chembiol.2023.05.001`), and `cftr-dms-tedman-2025` (DOI `10.7554/elife.107180`) remain `access: manual` / `claim_label: unverified` in `sources/catalog.yaml`; `load_mavedb_cc0_csv()` raises `PermissionError` until a human maintainer admits their license and terms.

## Split Hygiene & `PRED-0010` Outcome (`refuted`)

- **Positive controls (`7/7` PASS)**: Reproduces all 7 CFTR2 Variants of Unknown Significance in Table 1 of McDonald et al. (2024) (`R31L=0.17`, `V201M=0.36`, `A349V=0.26`, `S912L=0.12`, `D924N=0.83`, `M952T=0.85`, `Y1014C=0.37`), with `0` codon overlap against the 30 benchmark variants.
- **Domain- and position-disjoint split (`2026-11-21.v2`)**: `10` `train` variants in `{"Lasso_N_terminal", "TMD1", "R_domain"}` and `20` held-out `test` variants in unseen `{"NBD1", "TMD2", "NBD2_C_terminal"}` (`0` codon overlap, minimum cross-split codon gap `71 > 5` codons, `0` domain overlap), plus 4-fold Leave-One-Family-Out (`LOFO`) evaluation across `{"missense", "inframe_indel", "splice_regulatory", "class_i_ptc_frameshift"}`.
- **Recorded refutation (`PRED-0010`)**: Pre-registered responsiveness scoring `s_AM(x) = 1 - AM(x)` on the held-out 20-variant `TRIKAFTA` v23 `test` split achieved coverage `0.5500` (`11/20`), conditional missense ROC-AUC `0.5333` (`8/15`), and abstention-aware ROC-AUC `0.5101` (`101/198`), falling below the deterministic sequence-consequence class baseline (`1.0000` coverage, `0.6818` / `15/22` on `TRIKAFTA` v23; `0.8750` / `7/8` on `ALYFTREK` v9). Both modulator-responsive (`G551D=0.9897`, `N1303K=0.9907`, `S549N=0.9503`) and modulator-refractory (`L467P=0.9858`, `Y569D=0.9964`, `R560T=0.9846`) CF-causing missense variants receive high evolutionary pathogenicity scores, showing on open FDA SPL labels that sequence pathogenicity by itself does not predict pharmacological rescue.

## Verification Commands

```bash
python -m unittest discover -s benchmarks/theratype -v
python benchmarks/theratype/cftr_theratype_benchmark.py --check --falsify
python benchmarks/theratype/independent_audit.py --check
python ../cf-lab/.claude/skills/discovery-practice/scripts/check_claim_ledger.py --git benchmarks/theratype/ledger.json
```
