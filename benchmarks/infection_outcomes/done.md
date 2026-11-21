# Done file: CFF Infection Research Initiative NTM, phage, and non-sputum diagnostic trial outcomes benchmark (`T-0944` / `PRED-0009`)

Written and committed before running the benchmark or recording a claim ledger.

## Success test

The benchmark work is done when every condition below holds, each decided by an automated script:

1. `python -m unittest discover -s tools/trial_atlas -v` exits `0` from the repository root using only the Python 3.10+ standard library and local `tools/trial_atlas/` modules, with zero network calls and zero import dependency on sibling `cf-skills`.
2. `python tools/trial_atlas/infection_outcomes_benchmark.py --check` and `python tools/trial_atlas/infection_outcomes_benchmark.py --falsify` exit `0`, rejecting all five planted negative controls:
   - `OBSERVATIONAL` cohort leak (`NCT07312734` `SEND-CF` placed into the `940` PRISMA interventional cohort).
   - `X4`/`X5` non-PRISMA cohort leak (`NCT04294043` `ABATE` or `NCT05616221` non-CF bronchiectasis `AP-PA02` marked `included_in_prisma_940=True`).
   - Phase 1 primary-endpoint misattribution (claiming primary sputum bacterial `CFU/mL` density on safety-primary trials `NCT07275905` `AchromoPhage`, `NCT04596319` `AP-PA02`, `NCT05010577` `BX004-A`, or `NCT04294043` `ABATE`).
   - Trojan-Source bidi-split closing tag (`</untrusted\u200b_data_...>`), Cyrillic/Greek/fullwidth homoglyph `</untrusted_data_*>` breakout, spoofed IDN registry host, missing quarantine cage, or personal clinical eligibility overclaim.
   - `CC BY-NC-ND 4.0` derivative prose/table/translation leak on `10.1038/s41467-025-60598-4` or `10.1002/ppul.27475`, or citing `cff-research-priorities` (`unverified`) as numeric clinical proof.
3. On the fixed 20-study real-data Cystic Fibrosis infection and microbial-diagnostics held-out benchmark (`2026-11-21.v1`, disjoint from the 4 positive controls below), `infection_outcomes_benchmark.py --benchmark` achieves joint exact-match accuracy `>= 0.90` (`held_out_error <= 0.10` error fraction) across `(cohort_partition, pathogen_focus, endpoint_locus)` and beats the primary-outcome-only `PRISMA-940` baseline:
   - `NCT06998043`: `(prisma_940_interventional, bacteriophage_pseudomonas, primary_cfu_or_culture)`
   - `NCT05010577`: `(prisma_940_interventional, bacteriophage_pseudomonas, primary_safety_or_pk)`
   - `NCT04596319`: `(prisma_940_interventional, bacteriophage_pseudomonas, primary_safety_or_pk)`
   - `NCT04684641`: `(prisma_940_interventional, bacteriophage_pseudomonas, primary_cfu_or_culture)`
   - `NCT07275905`: `(prisma_940_interventional, bacteriophage_achromobacter, primary_safety_secondary_qpcr)`
   - `NCT05616221`: `(raw_1053_x5_excluded, non_cf_negation_control, excluded_non_cf)`
   - `NCT04294043`: `(raw_1053_x4_excluded, ntm_antimycobacterial, primary_safety_secondary_ntm_conversion)`
   - `NCT03597347`: `(prisma_940_interventional, ntm_antimycobacterial, primary_cfu_or_culture)`
   - `NCT05101915`: `(prisma_940_interventional, ntm_antimycobacterial, primary_cfu_or_culture)`
   - `NCT04685720`: `(prisma_940_interventional, ntm_antimycobacterial, primary_safety_or_pk)`
   - `NCT02354859`: `(prisma_940_interventional, pseudomonas_antimicrobial, primary_clinical_fev1)`
   - `NCT00721071`: `(prisma_940_interventional, non_sputum_diagnostic_sampling, primary_cfu_or_culture)`
   - `NCT01594827`: `(prisma_940_interventional, mrsa_staphylococcus, primary_cfu_or_culture)`
   - `NCT01746095`: `(prisma_940_interventional, mrsa_staphylococcus, primary_cfu_or_culture)`
   - `NCT03181932`: `(prisma_940_interventional, mrsa_staphylococcus, primary_clinical_fev1)`
   - `NCT00391976`: `(prisma_940_interventional, pseudomonas_antimicrobial, primary_cfu_or_culture)`
   - `NCT02453789`: `(prisma_940_interventional, burkholderia_gram_negative, primary_cfu_or_culture)`
   - `NCT06654752`: `(prisma_940_interventional, antibiotic_stewardship, primary_clinical_fev1)`
   - `NCT04545515`: `(prisma_940_interventional, modulator_negative_control, non_infection_control)`
   - `NCT07312734`: `(external_observational, non_sputum_diagnostic_sampling, observational_diagnostic_concordance)`
4. `python benchmarks/infection_outcomes/independent_audit.py --check` (regex-free token-window parser + `data/care_team_trials_directory.json` cross-verification, zero imports from `tools/trial_atlas/`) independently verifies `1,053` `INTERVENTIONAL` snapshot records, `940` PRISMA-included interventional studies, and `9` primary microbiology-outcome trials among the `17` PRISMA-940 benchmark studies with tolerance `0`.
5. `python tools/claims/check_anti_hallucination.py`, `python tools/claims/retraction_citation_radar.py --check`, `python ledger/check_predictions.py`, `python tools/sources/check_repo_organization.py`, and `python ../cf-lab/.claude/skills/discovery-practice/scripts/check_claim_ledger.py --git benchmarks/infection_outcomes/ledger.json` all exit `0`.

## Positive control

Before scoring any of the 20 held-out benchmark studies above, the classifier must reproduce four non-overlapping registry controls from `tools/trial_atlas/snapshot.json` (`NCT02212587`, `NCT01093521`, `NCT00910351`, `NCT01515007` — none of which appear in the 20-study benchmark):

1. `NCT02212587` (`TIP` Phase 1 *Burkholderia cepacia* complex) maps to `cohort_partition="prisma_940_interventional"` and `endpoint_locus="primary_cfu_or_culture"` (`primary_cfu_ml_or_g_reported=True`).
2. `NCT01093521` (`IV Ganite` Phase 1 gallium nitrate PK/safety in CF) maps to `cohort_partition="prisma_940_interventional"` and `endpoint_locus="primary_safety_or_pk"` (`primary_cfu_ml_or_g_reported=False`).
3. `NCT00910351` (`Cipro Inhaler` in CF children, `conditions=["Pseudomonas Infection"]`) maps to `cohort_partition="raw_1053_x4_excluded"` (`included_in_prisma_940=False` via rule `X4`).
4. `NCT01515007` (`ORBIT-3` Ciprofloxacin in Non Cystic Fibrosis Bronchiectasis) maps to `cohort_partition="raw_1053_x5_excluded"` (`included_in_prisma_940=False` via negation rule `X5`) and `endpoint_locus="excluded_non_cf"`.

If any positive control fails, `infection_outcomes_benchmark.py` raises `InfectionBenchmarkError` and reports no benchmark score.

## Stopping rule

Fixed 20-study held-out benchmark set (`2026-11-21.v1`) and 4 disjoint positive controls fixed in advance in this file; single evaluation pass with zero threshold or label tuning after seeing scores, and zero synthetic study records inside `--benchmark`. Every stub error, baseline miss, or review correction is logged in `benchmarks/infection_outcomes/failures.jsonl` and retained permanently.
