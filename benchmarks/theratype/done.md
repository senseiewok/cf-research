# Done file: Open CFTR theratype benchmark — does AlphaMissense predict CFTR modulator response from open variant-level data? (`T-0925` / `PRED-0010`)

Written and committed before running the `2026-11-21.v2` benchmark or recording a claim ledger.

## Success test

The benchmark work is done when every condition below holds, each decided by an automated script:

1. `python -m unittest discover -s benchmarks/theratype -v` exits `0` from the repository root using only the Python 3.10+ standard library and local `benchmarks/theratype/` files, with zero network calls and zero patient or PHI data.
2. `python benchmarks/theratype/cftr_theratype_benchmark.py --check --falsify` exits `0`, enforcing:
   - Zero synthetic or unmeasured model scores: every AlphaMissense, REVEL, CADD Phred, ClinVar Variation ID, and US FDA Structured Product Label (SPL) modulator-responsiveness value (`TRIKAFTA` NDA 212273/217660 v23 Table 6, `ALYFTREK` NDA 218730 v9 Table 5, `KALYDECO` NDA 203188 v36 Table 3, `SYMDEKO` NDA 210491 v22 Table 6) traces to admitted `sources/catalog.yaml` entries (`myvariant`, `openfda-label`, `clinvar`, `uniprot-p13569`, `europe-pmc` `PMC10810519`), with `None` (`"not stated"`) where unmeasured (`evo2_zero_shot_score=None`, `alphagenome_avi_score=None`, and `alphamissense_score=None` / `revel_score=None` on all non-missense indels, splicing, and Class I nonsense/frameshift alleles).
   - Deterministic structural-domain split auditing (`2026-11-21.v2`): `train` (`10` variants: `P67L`, `G85E`, `R117H`, `Y122X`, `621+1G>T`, `G178R`, `L206W`, `R334W`, `R347P`, `2051delAA` in `{"Lasso_N_terminal", "TMD1", "R_domain"}`) and `test` (`20` held-out variants: `A455E`, `L467P`, `I507del`, `F508del`, `V520F`, `1585-1G>A`, `G542X`, `S549N`, `G551D`, `R553X`, `R560T`, `Y569D`, `2789+5G>A`, `R1066C`, `D1152H`, `R1162X`, `3849+10kbC>T`, `S1251N`, `W1282X`, `N1303K` in unseen domains `{"NBD1", "TMD2", "NBD2_C_terminal"}`) have zero codon overlap, minimum cross-split codon distance `> 5` codons, zero recomputed `P13569` domain overlap, and zero `_heldout` family-renaming bypasses, alongside 4-fold Leave-One-Family-Out (`LOFO`) cross-validation across `{"missense", "inframe_indel", "splice_regulatory", "class_i_ptc_frameshift"}`.
   - Pre-registered directional responsiveness scoring for `PRED-0010`: using `s_AM(x) = 1.0 - alphamissense_score(x)` against target `Y = trikafta_v23_responsive` in `{0, 1}` (forbidding post-hoc sign flipping or `max(AUC, 1 - AUC)`), compute held-out `test` coverage (`|S|/N`), conditional missense ROC-AUC (`AUC_cond`), and abstention-aware ROC-AUC (`AUC_abst = 0.5 + (|S_1| * |S_0| / (N_1 * N_0)) * (AUC_cond - 0.5)`, where unscored pairs count as `0.5` ties), compared against a deterministic sequence-consequence class baseline (`stop_gained`/`frameshift`/canonical `+/-1,2` splice -> `0.0`; non-canonical intronic splice or in-frame indel -> `0.5`; missense -> `0.75`, fixed from `train` class biology without inspecting `test`). `PRED-0010` tests whether AlphaMissense achieves `test_coverage >= 0.80`, `conditional_missense_auroc >= 0.75`, and `abstention_aware_auroc >= 0.75` and matches or beats the sequence-consequence baseline; if any threshold fails, `PRED-0010` is recorded as `"refuted"`.
   - Rejection of all six planted negative controls:
     1. `codon_leakage_detected`: duplicate codon across `train` and `test`.
     2. `domain_or_buffer_leakage_detected`: cross-split domain overlap or codon distance `<= 5`.
     3. `heldout_suffix_hack_rejected`: non-canonical `variant_family` label (such as `gating_missense_heldout`) attempting to evade family checks.
     4. `non_missense_alphamissense_rejected`: non-`None` `alphamissense_score` or `revel_score` on an indel (`F508del`, `I507del`), nonsense (`Ter`), frameshift, or `c.` splice variant, or non-`None` uncataloged `evo2_zero_shot_score` / `alphagenome_avi_score`.
     5. `unadmitted_manual_source_rejected`: attempting to auto-ingest `access: manual` / `claim_label: unverified` sources (`mavedb-dump-zenodo`, `cftr-dms-mckee-2023`, `cftr-dms-tedman-2025`) before maintainer admission.
     6. `unknown_model_rejected`: unknown model identifier passed to the evaluator.
3. `python benchmarks/theratype/independent_audit.py --check` (separate code with zero imports from `cftr_theratype_benchmark.py` or `tools/`, independently reading `theratype_open_annotations_2026-11-21.json` and recomputing `P13569` domain spans, minimum cross-split codon gap, `TRIKAFTA` v23 vs `ALYFTREK` v9 responsive counts, AlphaMissense coverage, and exact tie-aware Mann-Whitney U ROC-AUC) agrees with tolerance `0` on integer counts and `1e-6` on floating-point metrics.
4. `python -m unittest discover -s tools/discovery -p "test_cftr_variant_to_delivery_deducer.py" -v`, `python tools/claims/check_anti_hallucination.py`, `python tools/discovery/discovery_council.py --falsify-all`, `python ledger/check_predictions.py`, `python tools/sources/check_repo_organization.py`, and `python ../cf-lab/.claude/skills/discovery-practice/scripts/check_claim_ledger.py --git benchmarks/theratype/ledger.json` all exit `0`.

## Positive control

Before scoring any of the 30 `train` + `test` benchmark variants, the evaluator must reproduce all seven CFTR2 Variants of Unknown Significance (VUS) reported in Table 1 (`pone.0297560.t001`) of McDonald et al., *PLOS ONE* 2024 (`10.1371/journal.pone.0297560`, `PMC10810519`, CC BY 4.0) from canonical `ENST00000003084` / `P13569` `myvariant` `dbnsfp.alphamissense` records (`0` codon overlap with all 30 `train` and `test` variants):

1. `R31L` (`p.Arg31Leu`, codon `31`): `round(score, 2) == 0.17`, class `"benign"`.
2. `V201M` (`p.Val201Met`, codon `201`): `round(score, 2) == 0.36`, class `"ambiguous"`.
3. `A349V` (`p.Ala349Val`, codon `349`): `round(score, 2) == 0.26`, class `"benign"`.
4. `S912L` (`p.Ser912Leu`, codon `912`): `round(score, 2) == 0.12`, class `"benign"`.
5. `D924N` (`p.Asp924Asn`, codon `924`): `round(score, 2) == 0.83`, class `"pathogenic"`.
6. `M952T` (`p.Met952Thr`, codon `952`): `round(score, 2) == 0.85`, class `"pathogenic"`.
7. `Y1014C` (`p.Tyr1014Cys`, codon `1014`): `round(score, 2) == 0.37`, class `"ambiguous"`.

Additionally, exact Mann-Whitney / abstention-ROC unit controls must confirm `AUC_abst == 0.5` on an all-`None` predictor, `1.0` on a strict positive separator, and `0.0` on an inverted separator (verifying zero post-hoc sign-flipping). If any positive control fails, `TheratypeBenchmarkError` is raised and no benchmark score is reported.

## Stopping rule

Fixed 30-variant benchmark (`10` `train` in `Lasso_N_terminal`/`TMD1`/`R_domain`, `20` `test` in `NBD1`/`TMD2`/`NBD2_C_terminal`, split `2026-11-21.v2`) and 7 disjoint McDonald et al. 2024 Table-1 positive controls fixed in advance in this file; single evaluation pass with zero split, threshold, or sign tuning after seeing held-out `test` scores, zero synthetic variant scores, and fail-closed refusal of `access: manual` sources (`mavedb-dump-zenodo`, `cftr-dms-mckee-2023`, `cftr-dms-tedman-2025`) until human maintainer ratification. Every defect found in the superseded `2026-11-15.v1` stub and any review corrections are logged in `benchmarks/theratype/failures.jsonl` and retained permanently.
