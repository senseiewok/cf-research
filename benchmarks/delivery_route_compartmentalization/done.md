# Done file: Delivery-route biomarker compartmentalization benchmark across early-phase CFTR interventional trials (`PRED-0015`)

Pre-registered and committed on `feat/sentinel-research-sync` prior to reading `tools/trial_atlas/snapshot.json` outcome counts or authoring `ledger.json`.
Independent open-source research benchmark only (Tier T1/T2 ClinicalTrials.gov registry metadata); not clinical guidance and not an assessment of personal trial eligibility or treatment efficacy. ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement.

## Success test

The benchmark work is complete when every deterministic check below exits `0` on the locked protocol without post-hoc modification of phase, route, regimen, or outcome rules:

1. **Dataset provenance, flattened schema check, & PRISMA-940 cohort boundary (`AH-1`, `AH-2`):**
   - Source dataset: `tools/trial_atlas/snapshot.json` (`catalog.yaml` source `clinicaltrials-gov`, `dataTimestamp` `2026-10-09`, `1,053` raw `INTERVENTIONAL` studies and `940` PRISMA-included Cystic Fibrosis interventional studies after `scope.py` rules `X1`–`X5`: `72` `X4` condition-omission exclusions and `41` `X5` negated-CF exclusions).
   - Explicitly assert that every PRISMA-included study record carries flattened `primary_outcomes` (`list` of dicts with `measure`, `description`, `time_frame`), `phases`, `interventions`, `brief_title`, and `official_title`; do not read uncollected `secondary_outcomes` or `outcomesModule` keys.

2. **Single pre-registered primary estimand (early-phase non-healthy-volunteer CFTR-restoring trials, zero OR-gates):**
   - **Early-phase filter:** `phases` (uppercase list) is non-empty and `set(phases) <= {"EARLY_PHASE1", "PHASE1", "PHASE2"}` (admitting `EARLY_PHASE1`, `PHASE1`, `["PHASE1", "PHASE2"]`, and `PHASE2`; excluding `["PHASE2", "PHASE3"]`, `PHASE3`, `PHASE4`, `NA`, and empty `[]`). All-phase counts across `940` PRISMA studies are reported strictly as non-gated descriptive context.
   - **Healthy-volunteer PK exclusion:** Exclude studies whose title (`brief_title` + `official_title`) names a healthy-volunteer cohort (`healthy volunteer(s)`, `healthy subject(s)`, `healthy participant(s)`, `healthy adult(s)`, `healthy male(s)`, `healthy female(s)`) without also naming people/patients/subjects/participants/adults with Cystic Fibrosis.
   - **Investigational-route precedence & CFTR-restoring two-arm partition:**
     - **Arm 2 (`inhaled_cftr` — airway-restricted respiratory-mucosal CFTR restoration, evaluated first):** Studies testing an investigational airway/respiratory-mucosal CFTR-restoring nucleic-acid, viral-vector, non-viral gene-transfer, splice-switching oligonucleotide, or topical-aminoglycoside nonsense-readthrough regimen:
       - Named inhaled/respiratory CFTR clinical candidates (`VX-522`, `MRT5005`, `ARCT-032`, `RCT2100`/`RCT-2100`, `BMD003`, `LUNAR-CF`, `4D-710`, `SP-101`, `KB407`, `BI 3720931`/`Lenticlair`, `pGM169`/`GL67A`, `eluforsen`/`QR-010`, `SPL84`/`SPL84-23-1`), OR
       - Explicit CFTR gene/mRNA/cDNA/plasmid/AAV/adenoviral/lentiviral/liposomal/antisense transfer paired with respiratory-mucosal delivery (`inhaled`, `inhalation`, `nebulized`, `nebulizer`, `nebuliser`, `aerosol`, `aerosolized`, `pulmonary`, `endobronchial`, `bronchoscopic`, `intranasal`, `nasal`, `maxillary sinus`, `sinus`, `airway`), OR
       - Topical/nebulized/intranasal `gentamicin`/`amikacin` explicitly paired with nonsense/stop-codon/PTC readthrough (`nonsense`, `stop codon`, `premature stop`, `readthrough`, `W1282X`, `G542X`).
       - Investigational-route precedence: co-administered background systemic modulators in an inhaled genetic trial (e.g., `VX-522` + `ivacaftor` in `NCT05668741`) keep the trial in `inhaled_cftr` with `dual_route_modulator_co_admin = True`, tested separately in a dual-route exclusion sensitivity check.
     - **Arm 1 (`systemic_cftr` — systemic small-molecule CFTR modulators & systemic PTC readthrough):** Studies not in `inhaled_cftr` whose `interventions[]` (`DRUG`/`BIOLOGICAL`/`GENETIC`) include a canonical systemic CFTR modulator (`ivacaftor`/`VX-770`/`Kalydeco`, `deutivacaftor`/`VX-561`/`CTP-656`, `lumacaftor`/`VX-809`/`Orkambi`, `tezacaftor`/`VX-661`/`Symdeko`/`Symkevi`, `elexacaftor`/`VX-445`/`Trikafta`/`Kaftrio`, `vanzacaftor`/`VX-121`/`Alyftrek`, `VX-659`, `VX-440`, `VX-152`, `VX-581`, `galicaftor`/`ABBV-2222`/`GLPG2222`, `navocaftor`/`ABBV-3067`/`GLPG3067`, `GLPG1837`, `GLPG2451`, `GLPG2737`, `ABBV-191`, `ABBV-576`, `posenacaftor`/`PTI-801`, `nesolicaftor`/`PTI-428`, `dirocaftor`/`PTI-808`, `FDL169`, `icenticaftor`/`QBW251`, `cavosonstat`/`N9115`, `4-phenylbutyrate`/`sodium phenylbutyrate`/`Buphenyl`) or systemic PTC readthrough agent (`ataluren`/`PTC124`, `ELX-02`), excluding therapy-withdrawal/discontinuation studies (`SIMPLIFY`, `"discontinuation"`, `"stopping"`, `"withdrawal of hypertonic"`) and extrinsic probe DDI studies (`"drug-drug interaction"`, `"midazolam"`, `"oral contraceptive"`, `"itraconazole"`, `"rifampin"`).
     - **Strict non-CFTR exclusion:** Topical non-CFTR ion-channel modulators (`ETD001`, `ARO-ENaC`, `VX-371`/`P-1037`, `SPX-101`, `BI 1265162`, `AZD5634`, `QBW276`, `amiloride`, `denufosol`), mucolytics/osmotics (`hypertonic saline`, `dornase alfa`, `mannitol`), and bacteriophage/antibiotic anti-infectives (`BX004`, `AP-PA02`, `AchromoPhage`, `tobramycin`, `aztreonam`, `colistin`, `gallium`) are excluded (`non_cftr_or_other`).

3. **Outcome definition & statistical decision criteria (`PRED-0015`):**
   - **Primary sweat chloride (`has_primary_sweat_chloride`):** At least one `primary_outcomes[].measure` string matches sweat chloride (`sweat chloride`, `sweat cl`, `sweat [cl-]`, `chloride concentration in sweat`) outside diagnostic/eligibility/inclusion wording (`diagnosis`, `diagnostic`, `inclusion`, `eligibility`, `screening criteria`).
   - **Pre-registered support criteria (all five required simultaneously for `standing = "open"` / `outcome = "supported"`; otherwise recorded as `refuted`):**
     1. Minimum early-phase arm margins: `N_sys >= 8` and `N_inh >= 8`.
     2. Directional compartmentalization: `p_sys = k_sys / N_sys > p_inh = k_inh / N_inh`.
     3. One-sided Fisher hypergeometric `p_fisher_1s < 0.01` (`Fraction(1, 100)`).
     4. Both the unconditional Chan-Zhang / Miettinen-Nurminen profile-score 95% two-sided risk-difference CI lower bound (`alpha/2 = 0.025`) `> 0.0` AND the Newcombe-Wilson Method-10 95% two-sided risk-difference CI lower bound `> 0.0`.
     5. Dual-route sensitivity preservation: after excluding `dual_route_modulator_co_admin == True` studies from `inhaled_cftr`, `p_sys > p_inh_pure` and one-sided Fisher `p_fisher_pure_1s < 0.01`.

4. **Independent second-method audit & repository gates (`discovery-practice` Rule 2):**
   - `benchmarks/delivery_route_compartmentalization/independent_audit.py` (importing zero `re`, zero `tools/`, and zero benchmark helper modules) independently verifies `done_sha256`, PRISMA `940` count parity against `data/care_team_trials_directory.json`, token-tuple arm & outcome counts, multiplicative `Fraction` hypergeometric `p_fisher_1s` (`tolerance = 0`), and 50-digit `Decimal` Newcombe-Wilson risk-difference bounds (`tolerance <= 1e-6`).
   - Unit tests (`test_delivery_route_compartmentalization.py`), `--check`, `--falsify`, `independent_audit.py --check`, `check_claim_ledger.py --git`, `ledger/check_predictions.py`, `check_anti_hallucination.py`, `check_repo_organization.py`, `check_codeql_patterns.py`, and `check_claim_strength.py` all exit `0`.

## Positive control

Before computing or inspecting observed early-phase sweat-chloride proportions on `snapshot.json`, `delivery_route_compartmentalization.py` must verify:

1. **Hypergeometric & CI calibration:** On synthetic fixture `(N_sys=12, k_sys=9, N_inh=10, k_inh=0)`, rational `p_fisher_1s == Fraction(11, 24871)`, full hypergeometric support sums to `Fraction(1, 1)`, and both Chan-Zhang unconditional 95% CI lower bound and Newcombe-Wilson 95% CI lower bound exceed `0.35`.
2. **Known trial regimen-route controls:** `NCT05668741` (`VX-522` + `ivacaftor`) classifies as `inhaled_cftr` with `dual_route_modulator_co_admin=True`; `NCT05248230` (`4D-710` AAV) and `NCT06429176` (`SPL84` ASO) classify as `inhaled_cftr` with `dual_route_modulator_co_admin=False`; `NCT00846781` (`denufosol` P2Y2) classifies as `non_cftr_or_other`.
3. **Inhaled airway/safety primary-endpoint specificity:** At least `80%` of early-phase `inhaled_cftr` trials register an airway (`FEV1`, `LCI`, `NPD`, mucociliary) or safety/tolerability primary endpoint in `primary_outcomes[].measure`.

## Stopping rule

Single pre-registered evaluation pass on frozen `tools/trial_atlas/snapshot.json` (`dataTimestamp` `2026-10-09`). Zero post-hoc tuning of phase list, healthy-volunteer filter, drug/vector allowlists, outcome keywords, or `alpha = 0.01` / `95%` CI thresholds after reading `snapshot.json`. Any positive-control failure raises `CompartmentalizationBenchmarkError` before reporting outcomes; if any of the five pre-registered support criteria fails, the prediction is recorded honestly as `refuted` in `results_2026-11-24.json`, `ledger.json`, `failures.jsonl`, and `ledger/predictions.jsonl`.
