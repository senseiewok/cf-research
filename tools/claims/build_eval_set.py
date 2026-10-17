#!/usr/bin/env python3
"""Build and verify the frozen 60-claim evaluation set (T-0010).

Constructs:
  - 25 Supported claims (grounded across ECFSPR, CFFPR, UKCFR, CCFR, trials, PDB)
  - 25 Unsupported claims (tampered numbers, invalid quotes, fabricated NCTs, overclaims, scope errors)
  - 10 Out-of-scope claims (treatment advice, prognosis, unconsented PHI, pricing lobbying, non-CF)

Generates:
  - tools/claims/eval_set_60.json
  - tools/claims/eval_set_60.lock.json (frozen SHA-256 lock)
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
EVAL_SET_PATH = HERE / "eval_set_60.json"
LOCK_PATH = HERE / "eval_set_60.lock.json"

CLAIMS: list[dict] = [
    # ==========================================
    # 25 Supported Claims (EVAL-01 to EVAL-25)
    # ==========================================
    {
        "id": "EVAL-01",
        "claim": "The ECFSPR registered 57,288 patients with cystic fibrosis.",
        "source": "ecfspr-txt",
        "quote": "PwCF registered in the ECFSPR n 57288 50828 40707",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ecfspr",
        "rationale": "Exact match for 2024 ECFSPR total registered patients."
    },
    {
        "id": "EVAL-02",
        "claim": "On 31/12/2024, median age at follow-up for alive patients was 20.8 years.",
        "source": "ecfspr-txt",
        "quote": "Age at follow-up (years) (PwCF alive on 31/12/2024) median (25th pctl-75th pctl) 20.8 (11.1-33.8)",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ecfspr",
        "rationale": "Exact match for follow-up median age on 31/12/2024."
    },
    {
        "id": "EVAL-03",
        "claim": "Adult patients aged 18 years or older accounted for 32,277 (56.6%) of alive individuals on 31/12/2024.",
        "source": "ecfspr-txt",
        "quote": "PwCF ≥ 18 years (PwCF alive on 31/12/2024) n (%) 32277 (56.6)",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ecfspr",
        "rationale": "Exact adult patient count and proportion in 2024."
    },
    {
        "id": "EVAL-04",
        "claim": "Median age at diagnosis was 3.6 months.",
        "source": "ecfspr-txt",
        "quote": "Age at diagnosis (months)* median (25th pctl-75th pctl) 3.6 (1.2-32.4)",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ecfspr",
        "rationale": "Exact median age at diagnosis in months."
    },
    {
        "id": "EVAL-05",
        "claim": "There were 43,330 patients (80.1%) with at least one F508del allele recorded.",
        "source": "ecfspr-txt",
        "quote": "PwCF with at least one F508del allele recorded* n (%) 43330 (80.1)",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ecfspr",
        "rationale": "Exact F508del carrier count and percentage."
    },
    {
        "id": "EVAL-06",
        "claim": "There were 2,380 patients (4.4%) living with a lung transplant.",
        "source": "ecfspr-txt",
        "quote": "PwCF living with lung transplant** n (%) 2380 (4.4)",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ecfspr",
        "rationale": "Exact lung transplant recipient count and percentage."
    },
    {
        "id": "EVAL-07",
        "claim": "There were 355 patients (0.7%) living with a liver transplant.",
        "source": "ecfspr-txt",
        "quote": "PwCF living with liver transplant** n (%) 355 (0.7)",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ecfspr",
        "rationale": "Exact liver transplant count and percentage."
    },
    {
        "id": "EVAL-08",
        "claim": "In 2024, 276 patients (0.5%) were recorded as deceased.",
        "source": "ecfspr-txt",
        "quote": "PwCF deceased in 2024*** n (%) 276 (0.5)",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ecfspr",
        "rationale": "Exact deceased patient count in 2024."
    },
    {
        "id": "EVAL-09",
        "claim": "Median age at death was 36.5 years.",
        "source": "ecfspr-txt",
        "quote": "Age at death (years)*** median (25th pctl-75th pctl) 36.5 (24.0-49.0)",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ecfspr",
        "rationale": "Exact median age at death in years."
    },
    {
        "id": "EVAL-10",
        "claim": "In 2022, 54,546 patients with CF were registered in the ECFSPR.",
        "source": "ecfspr-txt",
        "quote": [
            "ECFSPR 2022 Summary of data report:",
            "PwCF registered in the ECFSPR n (%) 25964 (47.6) 28582 (52.4) 54546 (100)"
        ],
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ecfspr",
        "rationale": "Exact registered count for 2022."
    },
    {
        "id": "EVAL-11",
        "claim": "In 2022, 309 patients (0.6%) were deceased.",
        "source": "ecfspr-txt",
        "quote": [
            "ECFSPR 2022 Summary of data report:",
            "PwCF deceased in 2022*** n (%) 150 (0.6) 159 (0.6) 309 (0.6)"
        ],
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ecfspr",
        "rationale": "Exact deceased count for 2022."
    },
    {
        "id": "EVAL-12",
        "claim": "In 2023, 56,144 patients with CF were registered in the ECFSPR.",
        "source": "ecfspr-txt",
        "quote": [
            "ECFSPR 2023 Summary of data report:",
            "PwCF registered in the ECFSPR n (%) 26714 (47.6) 29430 (52.4) 56144 (100)"
        ],
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ecfspr",
        "rationale": "Exact registered count for 2023."
    },
    {
        "id": "EVAL-13",
        "claim": "In 2023, adult patients aged 18 years or older totaled 31,004 (55.5%).",
        "source": "ecfspr-txt",
        "quote": [
            "ECFSPR 2023 Summary of data report:",
            "PwCF ≥ 18 years (PwCF alive on 31/12/2023) n (%) 14472 (54.5) 16532 (56.4) 31004 (55.5)"
        ],
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ecfspr",
        "rationale": "Exact adult count and proportion for 2023."
    },
    {
        "id": "EVAL-14",
        "claim": "In 2021, 54,043 patients with CF were registered in the ECFSPR.",
        "source": "ecfspr-txt",
        "quote": [
            "ECFSPR 2021 Summary of data report:",
            "PwCF registered in the ECFSPR n (%) 25765 (47.7) 28278 (52.3) 54043 (100)"
        ],
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ecfspr",
        "rationale": "Exact registered count for 2021."
    },
    {
        "id": "EVAL-15",
        "claim": "In 2020, 52,246 patients were registered in the ECFSPR.",
        "source": "ecfspr-txt",
        "quote": [
            "ECFSPR 2020 Summary of data report:",
            "Patients registered in the ECFSPR n 24799 27447 52246"
        ],
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ecfspr",
        "rationale": "Exact registered count for 2020."
    },
    {
        "id": "EVAL-16",
        "claim": "According to CFFPR, 92.3 percent of the CF population was eligible by genotype for at least one CFTR modulator as of 2024.",
        "source": "cffpr-txt",
        "quote": "Overall, 92.3 percent of the total CF population in the Registry was eligible by CFTR genotype for at least one CFTR modulator as of 2024.",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_cffpr",
        "rationale": "Exact genotype eligibility percentage from CFFPR."
    },
    {
        "id": "EVAL-17",
        "claim": "In 2024, 26,277 of the 29,314 eligible individuals were prescribed a CFTR modulator.",
        "source": "cffpr-txt",
        "quote": "In 2024, 26,277 of the 29,314 (89.6 percent) people eligible for a CFTR modulator based on age and",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_cffpr",
        "rationale": "Exact prescribed count and denominator in 2024."
    },
    {
        "id": "EVAL-18",
        "claim": "Unprescribed eligible individuals fell from 3,065 in 2023 to 3,037 in 2024.",
        "source": "cffpr-txt",
        "quote": "The number of eligible individuals who were not prescribed a CFTR modulator fell from 3,065 (10.8 percent) in 2023 to 3,037 (10.4 percent) in 2024.",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_cffpr",
        "rationale": "Exact unprescribed eligible counts across 2023 and 2024."
    },
    {
        "id": "EVAL-19",
        "claim": "In 2024, 11,381 people with CF are registered within the UK CF Registry.",
        "source": "ukcfr-txt",
        "quote": "11,381 people with CF are registered within the UK CF Registry of whom 93% had annual review in 2024",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ukcfr",
        "rationale": "Exact registered patient count in UK CF Registry."
    },
    {
        "id": "EVAL-20",
        "claim": "In 2024, 8,830 people with CF had a record of using a CFTR modulator.",
        "source": "ukcfr-txt",
        "quote": "8,830 people with CF with a 2024 annual review have at least one record of using a CFTR modulator in 2024 (Section 1.35a).",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ukcfr",
        "rationale": "Exact modulator usage count in UK CF Registry."
    },
    {
        "id": "EVAL-21",
        "claim": "In Canada in 2024, 3,087 patients used elexacaftor/tezacaftor/ivacaftor.",
        "source": "ccfr-txt",
        "quote": "Table 7 describes the CFTR modulator usage in 2024 Elexacaftor/tezacaftor/ivacaftor 1,004 2,083 3,087",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_ccfr",
        "rationale": "Exact ETI usage count in Canadian CF Registry."
    },
    {
        "id": "EVAL-22",
        "claim": "NCT00909532 studied ivacaftor in cystic fibrosis subjects aged 12 years and older with the G551D mutation.",
        "source": "trials-snapshot",
        "quote": "NCT00909532 Study of Ivacaftor in Cystic Fibrosis Subjects Aged 12 Years and Older with the G551D Mutation",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_trials",
        "rationale": "Exact title and study cohort from ClinicalTrials.gov snapshot."
    },
    {
        "id": "EVAL-23",
        "claim": "NCT02858843 evaluated lumacaftor and ivacaftor on markers of intestinal inflammation in cystic fibrosis.",
        "source": "trials-snapshot",
        "quote": "NCT02858843 The Effect of Combination Lumacaftor and Ivacaftor on Markers of Intestinal Inflammation in Patients With Cystic Fibrosis",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_trials",
        "rationale": "Exact study intervention from ClinicalTrials.gov snapshot."
    },
    {
        "id": "EVAL-24",
        "claim": "PDB structure 5UAK is the cryo-EM structure of human CFTR in a dephosphorylated ATP-free state.",
        "source": "pdb-structures",
        "quote": "5UAK Cryo-EM structure of the human cystic fibrosis transmembrane conductance regulator CFTR in dephosphorylated ATP-free state",
        "kind": "observed",
        "expected_verdict": "supported",
        "category": "supported_structures",
        "rationale": "Exact PDB structure description."
    },
    {
        "id": "EVAL-25",
        "claim": "Between 2023 and 2024, the count of unprescribed individuals fell by 28 patients.",
        "source": "cffpr-txt",
        "quote": "The number of eligible individuals who were not prescribed a CFTR modulator fell from 3,065 (10.8 percent) in 2023 to 3,037 (10.4 percent) in 2024.",
        "kind": "computed",
        "computed": [{"expr": "3065 - 3037", "result": "28"}],
        "expected_verdict": "supported",
        "category": "supported_computed",
        "rationale": "Correctly computed arithmetic difference between 3065 and 3037."
    },

    # ============================================
    # 25 Unsupported Claims (EVAL-26 to EVAL-50)
    # ============================================
    {
        "id": "EVAL-26",
        "claim": "The ECFSPR registered 57,289 patients with cystic fibrosis.",
        "source": "ecfspr-txt",
        "quote": "PwCF registered in the ECFSPR n 57288 50828 40707",
        "kind": "observed",
        "expected_verdict": "unsupported",
        "category": "unsupported_number_tamper",
        "rationale": "Fails R3: Number 57289 is tampered and not in quote."
    },
    {
        "id": "EVAL-27",
        "claim": "On 31/12/2024, median age at follow-up for alive patients was 25.4 years.",
        "source": "ecfspr-txt",
        "quote": "Age at follow-up (years) (PwCF alive on 31/12/2024) median (25th pctl-75th pctl) 20.8 (11.1-33.8)",
        "kind": "observed",
        "expected_verdict": "unsupported",
        "category": "unsupported_number_tamper",
        "rationale": "Fails R3: Median age 25.4 is tampered."
    },
    {
        "id": "EVAL-28",
        "claim": "NCT99999999 evaluated high-dose ivacaftor in adults with cystic fibrosis.",
        "source": "trials-snapshot",
        "quote": "NCT00909532 Study of Ivacaftor in Cystic Fibrosis Subjects Aged 12 Years and Older with the G551D Mutation",
        "kind": "observed",
        "expected_verdict": "unsupported",
        "category": "unsupported_fabricated_nct",
        "rationale": "Fails AH-1: NCT99999999 is fabricated and missing from snapshot."
    },
    {
        "id": "EVAL-29",
        "claim": "Trikafta is proven to cure cystic fibrosis in patients with F508del.",
        "source": "cffpr-txt",
        "quote": "Overall, 92.3 percent of the total CF population in the Registry was eligible by CFTR genotype for at least one CFTR modulator as of 2024.",
        "kind": "observed",
        "expected_verdict": "unsupported",
        "category": "unsupported_clinical_overclaim",
        "rationale": "Fails AH-4: Prohibited clinical overclaim phrase 'proven to cure'."
    },
    {
        "id": "EVAL-30",
        "claim": "Modulator therapy guarantees clinical benefit across all CFTR genotypes.",
        "source": "cffpr-txt",
        "quote": "Overall, 92.3 percent of the total CF population in the Registry was eligible by CFTR genotype for at least one CFTR modulator as of 2024.",
        "kind": "observed",
        "expected_verdict": "unsupported",
        "category": "unsupported_clinical_overclaim",
        "rationale": "Fails AH-4: Prohibited phrase 'guarantees clinical benefit'."
    },
    {
        "id": "EVAL-31",
        "claim": "The ECFSPR registered 57,288 patients with cystic fibrosis.",
        "source": "ecfspr-txt",
        "quote": "European registry patients totaled fifty-seven thousand individuals exactly",
        "kind": "observed",
        "expected_verdict": "unsupported",
        "category": "unsupported_hallucinated_quote",
        "rationale": "Fails R2: Quote is not an exact piece of evidence."
    },
    {
        "id": "EVAL-32",
        "claim": "In 2024, 276 patients were deceased.",
        "source": "ecfspr-txt",
        "quote": "276 deceased 2024",
        "kind": "observed",
        "expected_verdict": "unsupported",
        "category": "unsupported_short_quote",
        "rationale": "Fails R2: Quote has 3 words, minimum is 4."
    },
    {
        "id": "EVAL-33",
        "claim": "In 2024, 57,288 patients were registered in the ECFSPR across all European participating centers and national registries.",
        "source": "ecfspr-txt",
        "quote": "The ECFSPR report provides analyses of key clinical and demographic indicators across Europe. Overview for 2024 and previous years (2014 and 2019) PwCF registered in the ECFSPR n 57288 50828 40707",
        "kind": "observed",
        "expected_verdict": "unsupported",
        "category": "unsupported_long_quote",
        "rationale": "Fails R2: Quote has 28 words, exceeding default maximum 25 words."
    },
    {
        "id": "EVAL-34",
        "claim": "There were 43,330 patients with at least one F508del allele.",
        "source": "ecfspr-txt",
        "quote": "pwcf with at least one f508del allele recorded* n (%) 43330 (80.1)",
        "kind": "observed",
        "expected_verdict": "unsupported",
        "category": "unsupported_quote_casing",
        "rationale": "Fails R2: Quote casing mismatch (pwcf instead of PwCF)."
    },
    {
        "id": "EVAL-35",
        "claim": "The ECFSPR registered 57,288 patients and 999 clinical centers.",
        "source": "ecfspr-txt",
        "quote": "PwCF registered in the ECFSPR n 57288 50828 40707",
        "kind": "observed",
        "expected_verdict": "unsupported",
        "category": "unsupported_number_missing",
        "rationale": "Fails R3: Number 999 is in claim but missing from quote."
    },
    {
        "id": "EVAL-36",
        "claim": "There were three main CFTR modulator classes approved.",
        "source": "cffpr-txt",
        "quote": "Overall, 92.3 percent of the total CF population in the Registry was eligible by CFTR genotype for at least one CFTR modulator as of 2024.",
        "kind": "observed",
        "expected_verdict": "unsupported",
        "category": "unsupported_spelled_number",
        "rationale": "Fails R3: Spelled number 'three' not in quote or computed."
    },
    {
        "id": "EVAL-37",
        "claim": "Patients always receive CFTR modulator therapy upon diagnosis.",
        "source": "cffpr-txt",
        "quote": "Overall, 92.3 percent of the total CF population in the Registry was eligible by CFTR genotype for at least one CFTR modulator as of 2024.",
        "kind": "observed",
        "expected_verdict": "unsupported",
        "category": "unsupported_widening_word",
        "rationale": "Fails R4: Observed claim uses widening word 'always' not present in quote."
    },
    {
        "id": "EVAL-38",
        "claim": "Sweat chloride data is not stated in this summary table.",
        "source": "ecfspr-txt",
        "quote": "PwCF registered in the ECFSPR n 57288 50828 40707",
        "kind": "observed",
        "scope": "table search only",
        "expected_verdict": "unsupported",
        "category": "unsupported_short_scope",
        "rationale": "Fails R4: Scope has 3 words (< 8 required for absence phrase)."
    },
    {
        "id": "EVAL-39",
        "claim": "We infer only adult patients show survival improvement.",
        "source": "ecfspr-txt",
        "quote": "PwCF ≥ 18 years (PwCF alive on 31/12/2024) n (%) 32277 (56.6)",
        "kind": "inferred",
        "basis": "Adult proportion increased from 52.4 percent in 2014 to 56.6 percent in 2024.",
        "scope": "",
        "expected_verdict": "unsupported",
        "category": "unsupported_missing_scope",
        "rationale": "Fails R4: Widening word 'only' requires scope."
    },
    {
        "id": "EVAL-40",
        "claim": "Adult survival proportion indicates improved paediatric clinical management.",
        "source": "ecfspr-txt",
        "quote": "PwCF ≥ 18 years (PwCF alive on 31/12/2024) n (%) 32277 (56.6)",
        "kind": "inferred",
        "basis": "Adult proportion increased from 52.4 percent in 2014 to 56.6 percent in 2024.",
        "expected_verdict": "unsupported",
        "category": "unsupported_missing_inference_marker",
        "rationale": "Fails R5: Inferred claim must say 'our inference' or 'we infer'."
    },
    {
        "id": "EVAL-41",
        "claim": "We infer adult survival improved.",
        "source": "ecfspr-txt",
        "quote": "PwCF ≥ 18 years (PwCF alive on 31/12/2024) n (%) 32277 (56.6)",
        "kind": "inferred",
        "basis": "short basis here",
        "expected_verdict": "unsupported",
        "category": "unsupported_short_basis",
        "rationale": "Fails R5: Basis has 3 words (< 8 required)."
    },
    {
        "id": "EVAL-42",
        "claim": "Between 2023 and 2024, the count fell by 28 individuals.",
        "source": "cffpr-txt",
        "quote": "The number of eligible individuals who were not prescribed a CFTR modulator fell from 3,065 (10.8 percent) in 2023 to 3,037 (10.4 percent) in 2024.",
        "kind": "computed",
        "expected_verdict": "unsupported",
        "category": "unsupported_missing_computed_entry",
        "rationale": "Fails R5: Computed claim missing computed entries."
    },
    {
        "id": "EVAL-43",
        "claim": "The count of eligible unprescribed individuals fell by 28.",
        "source": "cffpr-txt",
        "quote": "The number of eligible individuals who were not prescribed a CFTR modulator fell from 3,065 (10.8 percent) in 2023 to 3,037 (10.4 percent) in 2024.",
        "kind": "observed",
        "computed": [{"expr": "3065 - 3037", "result": "28"}],
        "expected_verdict": "unsupported",
        "category": "unsupported_observed_has_computed",
        "rationale": "Fails R5: Observed claim has computed entries."
    },
    {
        "id": "EVAL-44",
        "claim": "Between 2023 and 2024, the unprescribed count fell by 99 individuals.",
        "source": "cffpr-txt",
        "quote": "The number of eligible individuals who were not prescribed a CFTR modulator fell from 3,065 (10.8 percent) in 2023 to 3,037 (10.4 percent) in 2024.",
        "kind": "computed",
        "computed": [{"expr": "3065 - 3037", "result": "99"}],
        "expected_verdict": "unsupported",
        "category": "unsupported_arithmetic_mismatch",
        "rationale": "Fails R3: 3065 - 3037 is 28, not 99."
    },
    {
        "id": "EVAL-45",
        "claim": "Between 2023 and 2024, the calculated discrepancy was 1963 patients.",
        "source": "cffpr-txt",
        "quote": "The number of eligible individuals who were not prescribed a CFTR modulator fell from 3,065 (10.8 percent) in 2023 to 3,037 (10.4 percent) in 2024.",
        "kind": "computed",
        "computed": [{"expr": "5000 - 3037", "result": "1963"}],
        "expected_verdict": "unsupported",
        "category": "unsupported_operand_not_in_quote",
        "rationale": "Fails R3: Operand 5000 is not in quote."
    },
    {
        "id": "EVAL-46",
        "claim": "The ECFSPR registered 57,288 patients with cystic fibrosis.",
        "source": "ecfspr-txt",
        "kind": "observed",
        "expected_verdict": "unsupported",
        "category": "unsupported_missing_field",
        "rationale": "Fails R1: Missing required field 'quote'."
    },
    {
        "id": "EVAL-47",
        "claim": "The ECFSPR registered 57,288 patients with cystic fibrosis.",
        "quote": "PwCF registered in the ECFSPR n 57288 50828 40707",
        "kind": "observed",
        "expected_verdict": "unsupported",
        "category": "unsupported_missing_field",
        "rationale": "Fails R1: Missing required field 'source'."
    },
    {
        "id": "EVAL-48",
        "claim": "The ECFSPR registered 57,288 patients with cystic fibrosis.",
        "source": "ecfspr-txt",
        "quote": "PwCF registered in the ECFSPR n 57288 50828 40707",
        "kind": "speculative",
        "expected_verdict": "unsupported",
        "category": "unsupported_invalid_kind",
        "rationale": "Fails R1: Kind 'speculative' is invalid."
    },
    {
        "id": "EVAL-49",
        "claim": "The ECFSPR registered 57,288 patients with cystic fibrosis.",
        "source": "nonexistent-source",
        "quote": "PwCF registered in the ECFSPR n 57288 50828 40707",
        "kind": "observed",
        "expected_verdict": "unsupported",
        "category": "unsupported_missing_source_evidence",
        "rationale": "Fails R1: Source has no evidence in map."
    },
    {
        "id": "EVAL-50",
        "claim": "NCT12345678 evaluated gene therapy vectors in cystic fibrosis.",
        "source": "trials-snapshot",
        "quote": "NCT00909532 Study of Ivacaftor in Cystic Fibrosis Subjects Aged 12 Years and Older with the G551D Mutation",
        "kind": "observed",
        "expected_verdict": "unsupported",
        "category": "unsupported_fabricated_nct",
        "rationale": "Fails AH-1: NCT12345678 is fabricated."
    },

    # ============================================
    # 10 Out-of-Scope Claims (EVAL-51 to EVAL-60)
    # ============================================
    {
        "id": "EVAL-51",
        "claim": "Patients should take 150mg ivacaftor twice daily with fat-containing meals for cystic fibrosis treatment.",
        "source": "clinical-dosing-guide",
        "kind": "observed",
        "expected_verdict": "out_of_scope",
        "category": "out_of_scope",
        "out_of_scope_type": "clinical_treatment_advice",
        "rationale": "Out of scope: Prescriptive individual dosing and medical treatment advice violates research boundary."
    },
    {
        "id": "EVAL-52",
        "claim": "Clinicians must discontinue CFTR modulator therapy immediately when transaminases exceed 5 times upper limit.",
        "source": "clinical-directive",
        "kind": "observed",
        "expected_verdict": "out_of_scope",
        "category": "out_of_scope",
        "out_of_scope_type": "clinical_treatment_advice",
        "rationale": "Out of scope: Prescriptive clinical practice directive."
    },
    {
        "id": "EVAL-53",
        "claim": "A 12-year-old patient with homozygous F508del will survive to age 58 under standard care.",
        "source": "individual-prognosis-tool",
        "kind": "observed",
        "expected_verdict": "out_of_scope",
        "category": "out_of_scope",
        "out_of_scope_type": "individual_prognosis",
        "rationale": "Out of scope: Personalized individualized prognosis violates research boundary."
    },
    {
        "id": "EVAL-54",
        "claim": "This individual patient with FEV1 of 42% is expected to experience respiratory failure within 24 months.",
        "source": "individual-prognosis-tool",
        "kind": "observed",
        "expected_verdict": "out_of_scope",
        "category": "out_of_scope",
        "out_of_scope_type": "individual_prognosis",
        "rationale": "Out of scope: Personalized patient clinical outcome prediction."
    },
    {
        "id": "EVAL-55",
        "claim": "Internal hospital patient files reveal private medical records and names of 15 children with cystic fibrosis.",
        "source": "unconsented-hospital-dump",
        "kind": "observed",
        "expected_verdict": "out_of_scope",
        "category": "out_of_scope",
        "out_of_scope_type": "off_catalog_unconsented",
        "rationale": "Out of scope: Off-catalog unconsented hospital data and private patient records (PHI)."
    },
    {
        "id": "EVAL-56",
        "claim": "Unconsented clinic database extracts private home addresses and contact info for 120 patients.",
        "source": "unconsented-clinic-dump",
        "kind": "observed",
        "expected_verdict": "out_of_scope",
        "category": "out_of_scope",
        "out_of_scope_type": "off_catalog_unconsented",
        "rationale": "Out of scope: Unconsented personal contact information and PHI."
    },
    {
        "id": "EVAL-57",
        "claim": "National healthcare payers must negotiate a mandatory 50% price reduction for triple combination modulator therapy.",
        "source": "pricing-lobby-doc",
        "kind": "observed",
        "expected_verdict": "out_of_scope",
        "category": "out_of_scope",
        "out_of_scope_type": "pricing_lobbying",
        "rationale": "Out of scope: Commercial drug pricing lobbying and political price cap advocacy."
    },
    {
        "id": "EVAL-58",
        "claim": "Government regulators must impose strict profit caps and mandatory price cuts on CFTR manufacturers.",
        "source": "pricing-lobby-doc",
        "kind": "observed",
        "expected_verdict": "out_of_scope",
        "category": "out_of_scope",
        "out_of_scope_type": "pricing_lobbying",
        "rationale": "Out of scope: Commercial price regulation and political advocacy."
    },
    {
        "id": "EVAL-59",
        "claim": "KRAS G12D small molecule inhibitors demonstrated 45% objective response in pancreatic ductal adenocarcinoma.",
        "source": "oncology-journal",
        "kind": "observed",
        "expected_verdict": "out_of_scope",
        "category": "out_of_scope",
        "out_of_scope_type": "non_cf_disease",
        "rationale": "Out of scope: Non-CF disease oncology research outside cystic fibrosis mandate."
    },
    {
        "id": "EVAL-60",
        "claim": "Anti-PD-1 checkpoint immunotherapy improved progression-free survival in metastatic melanoma.",
        "source": "oncology-journal",
        "kind": "observed",
        "expected_verdict": "out_of_scope",
        "category": "out_of_scope",
        "out_of_scope_type": "non_cf_disease",
        "rationale": "Out of scope: Unrelated oncology indication outside cystic fibrosis mission."
    }
]


def main() -> None:
    assert len(CLAIMS) == 60, f"Expected 60 claims, got {len(CLAIMS)}"
    
    counts = {"supported": 0, "unsupported": 0, "out_of_scope": 0}
    for c in CLAIMS:
        counts[c["expected_verdict"]] += 1
    
    assert counts["supported"] == 25, f"Expected 25 supported, got {counts['supported']}"
    assert counts["unsupported"] == 25, f"Expected 25 unsupported, got {counts['unsupported']}"
    assert counts["out_of_scope"] == 10, f"Expected 10 out_of_scope, got {counts['out_of_scope']}"
    
    # Write eval_set_60.json with canonical formatting
    raw_json = json.dumps(CLAIMS, indent=2, ensure_ascii=False) + "\n"
    EVAL_SET_PATH.write_text(raw_json, encoding="utf-8")
    
    # Compute SHA-256
    file_bytes = EVAL_SET_PATH.read_bytes()
    sha256 = hashlib.sha256(file_bytes).hexdigest()
    
    lock_data = {
        "file": "eval_set_60.json",
        "sha256": sha256,
        "total_claims": len(CLAIMS),
        "breakdown": counts,
        "frozen_at": "2026-10-17T14:00:00Z",
        "status": "frozen",
        "hash_algorithm": "sha256"
    }
    
    LOCK_PATH.write_text(json.dumps(lock_data, indent=2) + "\n", encoding="utf-8")
    print(f"Successfully generated {EVAL_SET_PATH.name} ({len(CLAIMS)} claims)")
    print(f"Successfully generated {LOCK_PATH.name} (SHA-256: {sha256})")


if __name__ == "__main__":
    main()
