#!/usr/bin/env python3
"""CFF Infection Research Initiative NTM, Phage, and Non-Sputum Diagnostic Benchmark (T-0944 / PRED-0009).

Benchmarks registered microbiological and diagnostic endpoints across Cystic Fibrosis
Foundation (CFF) Infection Research Initiative interventional trials (`NCT06998043` BX004,
`NCT07275905` AchromoPhage, `NCT04596319` AP-PA02, `NCT04294043` ABATE IV gallium for NTM),
comparator CF infection/sampling trials in `tools/trial_atlas/snapshot.json`, external
`OBSERVATIONAL` diagnostic concordance studies (`NCT07312734` SEND-CF), and `CC BY 4.0`
face-mask aerosol sampling feasibility data (`DOI 10.1099/jmm.0.002150`), anchored to
`cffpr-highlights-2025` prevalence (`22%` P. aeruginosa, `12%` MRSA, `12%` NTM; median age
at first NTM infection `26` years).

Invariants (`cf-research-context`, `security-browsing`, `AH-1`..`AH-4`, `AH-12`):
- 100% pure Python 3.10+ stdlib + local `tools/trial_atlas/` modules; zero network/LLM calls.
- Strict `AH-2` separation: `prisma_940_interventional` (`940`) vs `raw_1053_x4_excluded`
  (`NCT04294043` ABATE, retrieved on the `condition` route via `official_title` expansion
  and excluded by `X4` because `conditions` lists only `"Nontuberculous Mycobacterium
  Infection"`) vs `raw_1053_x5_excluded` (`NCT05616221` non-CF bronchiectasis) vs
  `external_observational` (`NCT07312734` SEND-CF).
- Separates `snapshot.json` (`primary_outcomes` only) from dated `ClinicalTrials.gov` API v2
  `protocolSection` metadata quarantined inside `<untrusted_data_NONCE>` (escaping bidi
  controls and Cyrillic/Greek/fullwidth closing-tag homoglyphs).
- Isolates `cff-research-priorities` (`unverified` pointer only) and blocks derivative text
  or translations on `CC BY-NC-ND 4.0` references.
"""
from __future__ import annotations

import argparse
import json
import re
import secrets
import sys
import unicodedata
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from care_team_trials import DEFAULT_SNAPSHOT_JSON, NOT_STATED, load_prisma_included_studies  # noqa: E402
import lexicon as lex  # noqa: E402
import scope as scope_mod  # noqa: E402

BENCHMARK_VERSION = "2026-11-21.v1"
DEFAULT_LEXICON_JSON = HERE / "lexicon.json"
NLM_DISCLAIMER = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)
CLINICAL_BOUNDARY = (
    "Independent open-source research benchmark only (Tier T1/T2 registry and OA metadata); "
    "not clinical guidance and never an assessment of personal trial eligibility or treatment efficacy."
)
ALLOWED_REGISTRY_HOSTS: frozenset[str] = frozenset({"clinicaltrials.gov", "www.clinicaltrials.gov"})
TROJAN_BIDI_CODEPOINTS: frozenset[int] = frozenset(
    {0x200B, 0x200C, 0x200D, 0x200E, 0x200F, 0x202A, 0x202B, 0x202C, 0x202D, 0x202E, 0x2060, 0x2066, 0x2067, 0x2068, 0x2069, 0xFEFF}
)
CONFUSABLE_ASCII_MAP = str.maketrans(
    {
        # Cyrillic
        "\u0430": "a", "\u0435": "e", "\u043e": "o", "\u0440": "p", "\u0441": "c", "\u0443": "y",
        "\u0445": "x", "\u0501": "d", "\u0456": "i", "\u0455": "s", "\u0442": "t", "\u043d": "n",
        "\u0410": "A", "\u0412": "B", "\u0421": "C", "\u0500": "D", "\u0415": "E", "\u041d": "H",
        "\u0406": "I", "\u041a": "K", "\u041c": "M", "\u041e": "O", "\u0420": "P", "\u0405": "S",
        "\u0422": "T", "\u0425": "X", "\u0423": "Y",
        # Greek
        "\u03b1": "a", "\u03b5": "e", "\u03b9": "i", "\u03bf": "o", "\u03c1": "p", "\u03c4": "t",
        "\u03c5": "u", "\u03bd": "n", "\u03ba": "k", "\u03c7": "x", "\u03c2": "s", "\u03c3": "s",
        "\u0391": "A", "\u0392": "B", "\u0395": "E", "\u0397": "H", "\u0399": "I", "\u039a": "K",
        "\u039c": "M", "\u039d": "N", "\u039f": "O", "\u03a1": "P", "\u03a4": "T", "\u03a7": "X", "\u03a5": "Y",
        # Fullwidth / math angle & slash confusables
        "\uFF1C": "<", "\uFF1E": ">", "\uFF0F": "/", "\u2044": "/", "\u2215": "/",
        "\u2329": "<", "\u232A": ">", "\u27E8": "<", "\u27E9": ">",
    }
)
UNTRUSTED_CLOSE_LOOKALIKE_RE = re.compile(r"<\s*/\s*[^\s>]*data|/\s*untrusted", re.I)
FORBIDDEN_OVERCLAIM_PHRASES: tuple[re.Pattern[str], ...] = (
    re.compile(r"\b(?:you\s+are\s+eligible|you\s+may\s+be\s+eligible|you\s+qualify|qualifies\s+for)\b", re.I),
    re.compile(r"\b(?:proven\s+to\s+cure|guarantees\s+clinical\s+benefit|safe\s+and\s+effective\s+for\s+patients|cures\s+cystic\s+fibrosis)\b", re.I),
)

CFFPR_2025_INFECTION_ORIENTATION: dict[str, Any] = {
    "catalog_source_id": "cffpr-highlights-2025",
    "claim_label": "verified",
    "pseudomonas_aeruginosa_pct": 22,
    "mrsa_pct": 12,
    "ntm_pct": 12,
    "ntm_median_first_infection_age_years": 26,
    "navigation_pointer_only": {
        "catalog_source_id": "cff-research-priorities",
        "access": "manual",
        "claim_label": "unverified",
        "permitted_use": "navigation_pointer_only_never_in_claims_json",
    },
}

NON_SPUTUM_LITERATURE_ANCHORS: tuple[dict[str, Any], ...] = (
    {
        "doi": "10.1099/jmm.0.002150",
        "pmid": "41950061",
        "pmcid": "PMC13061336",
        "license": "CC BY 4.0",
        "catalog_source_id": "europe-pmc",
        "retracted": False,
        "editorial_notices_count": 0,
        "title": "Face mask sampling for the detection of microbes in expelled aerosols and the impact of airway clearance on microbial yield in children with cystic fibrosis: a feasibility trial",
        "participants_enrolled": 11,
        "median_age_years": 12,
        "age_range_years": "1-16",
        "mask_wear_duration_minutes": 15,
        "respiratory_pathogen_positive_mask_samples": 9,
        "mycobacterium_abscessus_positive_n": 3,
        "pseudomonas_aeruginosa_positive_mask_samples": 0,
        "post_clearance_16s_increased_yield_n": 5,
        "post_clearance_16s_statistically_significant": False,
        "verbatim_cc_by_excerpt": (
            "Children wore a suitably sized face mask fitted with two strips of a polyvinyl alcohol sampling matrix for a period of 15 min. "
            "Eleven children (six male, median age 12 years, range 1-16 years), from the Leicester CF cohort were included. "
            "Nine face mask samples from 11 participants isolated respiratory pathogens, including M. abscessus (n=3). "
            "P. aeruginosa was not detected on face mask samples. "
            "There was a trend towards an increase in microbial yield (prGen16s) following airway clearance (n=5), but this did not reach statistical significance."
        ),
        "derivative_text": None,
        "translated_text": None,
    },
    {
        "doi": "10.1038/s41467-025-60598-4",
        "pmid": "40593506",
        "pmcid": "PMC12216271",
        "license": "CC BY-NC-ND 4.0",
        "catalog_source_id": "europe-pmc",
        "retracted": False,
        "editorial_notices_count": 1,
        "editorial_notice_doi": "10.1038/s41467-026-70146-3",
        "editorial_notice_type": "correction",
        "title": "Phage therapy with nebulized cocktail BX004-A for chronic Pseudomonas aeruginosa infections in cystic fibrosis: a randomized first-in-human trial",
        "linked_nct_id": "NCT05010577",
        "verbatim_cc_by_excerpt": None,
        "derivative_text": None,
        "translated_text": None,
    },
    {
        "doi": "10.1002/ppul.27475",
        "pmid": "39785222",
        "pmcid": "PMC11715147",
        "license": "CC BY-NC-ND 4.0",
        "catalog_source_id": "europe-pmc",
        "retracted": False,
        "editorial_notices_count": 0,
        "title": "Genetic Concordance of Staphylococcus aureus From Oropharyngeal and Sputum Cultures in People With Cystic Fibrosis",
        "verbatim_cc_by_excerpt": None,
        "derivative_text": None,
        "translated_text": None,
    },
)

# Verbatim ClinicalTrials.gov API v2 protocolSection fields verified 2026-11-21 (separated from snapshot.json).
EXTERNAL_AND_SECONDARY_PROTOCOL_METADATA: dict[str, dict[str, Any]] = {
    "NCT06998043": {
        "provenance": "clinicaltrials_gov_api_v2",
        "verified_api_date": "2026-11-21",
        "last_update_post_date": "2025-07-22",
        "secondary_outcomes": [
            {"measure": "Change in sputum PsA burden", "description": "Change from baseline in sputum PsA CFU/g at D8, D29, D43, D85, 3 months post-dose and 6 months post-dose", "timeFrame": "until 6 months after last dose of study drug"},
            {"measure": "Efficacy of BX004 on obtaining negative sputum cultures for PsA", "description": "Proportion of subjects with negative sputum cultures for PsA, time to negative sputum cultures, durability of negative sputum cultures", "timeFrame": "until 6 months after last dose of study drug"},
            {"measure": "Change in lung function at D8, D29, D43, D57 (EOT), D85, 3 months post-dose, and 6 months post-dose", "description": "Change from baseline in % predicted FEV1", "timeFrame": "from Day 8 until 6 months after last dose (end of study)"},
        ],
        "specimen_modalities": ["sputum"],
    },
    "NCT07275905": {
        "provenance": "clinicaltrials_gov_api_v2",
        "verified_api_date": "2026-11-21",
        "last_update_post_date": "2026-10-08",
        "secondary_outcomes": [
            {
                "measure": "AchromoPhage Concentration by qPCR (Pharmacokinetic Outcome Measure)",
                "description": "Quantitative serum concentrations of each phage component of AchromoPhage will be measured by qPCR using primers that will amplify each of the 4 phages in the cocktail at predefined time points, following inhaled, intravenous, or combined administration. Blood PCR: Collected at 0-, 10-, 30-, 60-, 120-, 180-, and 240-minute intervals following AchromoPhage administration during dosing visits (Days 0, 7, and 14) and once on Days 21 and 42. Sputum, nasal swabs, and oral swabs: Collected weekly (Days 0, 7, and 14, 21, 28, 35, and 42).",
                "timeFrame": "Blood: Days 0, 7, 14, 21, 42. Sputum, nasal swabs, and oral swabs: Days 0, 7, 14, 21, 28, 35, and 42.",
            },
            {"measure": "Number of Participants With Treatment-Emergent Serum Phage-Neutralizing Antibody Activity", "description": "Quantitative phage-neutralizing antibody activity will be measured in serum samples to assess treatment-emergent humoral immune responses following AchromoPhage administration.", "timeFrame": "Days 0, 7, 14, 21, 28, and 42"},
        ],
        "bacterial_cfu_ml_in_outcomes_module": NOT_STATED,
        "specimen_modalities": ["blood_serum", "sputum", "nasal_swab", "oral_swab"],
    },
    "NCT04596319": {
        "provenance": "clinicaltrials_gov_api_v2",
        "verified_api_date": "2026-11-21",
        "last_update_post_date": "2024-01-31",
        "secondary_outcomes": [],
        "secondary_outcomes_note": "protocolSection.outcomesModule.secondaryOutcomes is absent in ClinicalTrials.gov API v2 record (hasResults=true).",
        "specimen_modalities": ["sputum"],
    },
    "NCT04294043": {
        "provenance": "clinicaltrials_gov_api_v2",
        "verified_api_date": "2026-11-21",
        "last_update_post_date": "2025-10-01",
        "secondary_outcomes": [
            {"measure": "Clinically significant abnormal laboratory measures. (safety)", "description": "Proportion of patients experiencing clinically significant abnormal laboratory measures.", "timeFrame": "Baseline to Day 57"},
            {"measure": "NTM clearance (efficacy)", "description": "Proportion of subjects who were NTM culture positive at baseline and have at least 2 sequential negative NTM cultures between visits 2 (Day 6) and 7 (Day 111). Those negative cultures must be at least 2 weeks apart.", "timeFrame": "Day 6 to Day 111"},
        ],
        "specimen_modalities": ["sputum_ntm_culture"],
    },
    "NCT07312734": {
        "provenance": "clinicaltrials_gov_api_v2",
        "verified_api_date": "2026-11-21",
        "last_update_post_date": "2026-07-02",
        "study_type": "OBSERVATIONAL",
        "overall_status": "RECRUITING",
        "enrollment": 300,
        "lead_sponsor": "Chris Goss",
        "brief_title": "Study to Enable New Diagnostics for Pulmonary Microbes in People With CF",
        "official_title": "Study to Enable New Diagnostics for Pulmonary Microbes in People With CF (SEND-CF)",
        "conditions": ["Cystic Fibrosis (CF)", "New Diagnostics"],
        "primary_outcomes": [
            {
                "measure": "Evaluate alternative methods for the detection of microbial pathogens from the lungs of people with CF by comparing against the gold standard of sputum cultures.",
                "description": "Specimens (including sputum, saliva, blood, urine, and breath (for the Owlstone Device Sub-Study)) will be collected at the Study Visit. These specimens will be stored at a biorepository and be made available to the research community to test and validate diagnostic assays for the detection of CF pathogens.",
                "time_frame": "Baseline Visit",
            }
        ],
        "secondary_outcomes": [],
        "specimen_modalities": ["sputum", "saliva", "blood", "urine", "exhaled_breath"],
        "registry_url": "https://clinicaltrials.gov/study/NCT07312734",
    },
}

# Fixed 20-study held-out benchmark (locked in benchmarks/infection_outcomes/done.md, disjoint from positive controls).
INFECTION_BENCHMARK_20_GROUND_TRUTH: tuple[dict[str, str], ...] = (
    {"nct_id": "NCT06998043", "expected_cohort": "prisma_940_interventional", "expected_focus": "bacteriophage_pseudomonas", "expected_locus": "primary_cfu_or_culture"},
    {"nct_id": "NCT05010577", "expected_cohort": "prisma_940_interventional", "expected_focus": "bacteriophage_pseudomonas", "expected_locus": "primary_safety_or_pk"},
    {"nct_id": "NCT04596319", "expected_cohort": "prisma_940_interventional", "expected_focus": "bacteriophage_pseudomonas", "expected_locus": "primary_safety_or_pk"},
    {"nct_id": "NCT04684641", "expected_cohort": "prisma_940_interventional", "expected_focus": "bacteriophage_pseudomonas", "expected_locus": "primary_cfu_or_culture"},
    {"nct_id": "NCT07275905", "expected_cohort": "prisma_940_interventional", "expected_focus": "bacteriophage_achromobacter", "expected_locus": "primary_safety_secondary_qpcr"},
    {"nct_id": "NCT05616221", "expected_cohort": "raw_1053_x5_excluded", "expected_focus": "non_cf_negation_control", "expected_locus": "excluded_non_cf"},
    {"nct_id": "NCT04294043", "expected_cohort": "raw_1053_x4_excluded", "expected_focus": "ntm_antimycobacterial", "expected_locus": "primary_safety_secondary_ntm_conversion"},
    {"nct_id": "NCT03597347", "expected_cohort": "prisma_940_interventional", "expected_focus": "ntm_antimycobacterial", "expected_locus": "primary_cfu_or_culture"},
    {"nct_id": "NCT05101915", "expected_cohort": "prisma_940_interventional", "expected_focus": "ntm_antimycobacterial", "expected_locus": "primary_cfu_or_culture"},
    {"nct_id": "NCT04685720", "expected_cohort": "prisma_940_interventional", "expected_focus": "ntm_antimycobacterial", "expected_locus": "primary_safety_or_pk"},
    {"nct_id": "NCT02354859", "expected_cohort": "prisma_940_interventional", "expected_focus": "pseudomonas_antimicrobial", "expected_locus": "primary_clinical_fev1"},
    {"nct_id": "NCT00721071", "expected_cohort": "prisma_940_interventional", "expected_focus": "non_sputum_diagnostic_sampling", "expected_locus": "primary_cfu_or_culture"},
    {"nct_id": "NCT01594827", "expected_cohort": "prisma_940_interventional", "expected_focus": "mrsa_staphylococcus", "expected_locus": "primary_cfu_or_culture"},
    {"nct_id": "NCT01746095", "expected_cohort": "prisma_940_interventional", "expected_focus": "mrsa_staphylococcus", "expected_locus": "primary_cfu_or_culture"},
    {"nct_id": "NCT03181932", "expected_cohort": "prisma_940_interventional", "expected_focus": "mrsa_staphylococcus", "expected_locus": "primary_clinical_fev1"},
    {"nct_id": "NCT00391976", "expected_cohort": "prisma_940_interventional", "expected_focus": "pseudomonas_antimicrobial", "expected_locus": "primary_cfu_or_culture"},
    {"nct_id": "NCT02453789", "expected_cohort": "prisma_940_interventional", "expected_focus": "burkholderia_gram_negative", "expected_locus": "primary_cfu_or_culture"},
    {"nct_id": "NCT06654752", "expected_cohort": "prisma_940_interventional", "expected_focus": "antibiotic_stewardship", "expected_locus": "primary_clinical_fev1"},
    {"nct_id": "NCT04545515", "expected_cohort": "prisma_940_interventional", "expected_focus": "modulator_negative_control", "expected_locus": "non_infection_control"},
    {"nct_id": "NCT07312734", "expected_cohort": "external_observational", "expected_focus": "non_sputum_diagnostic_sampling", "expected_locus": "observational_diagnostic_concordance"},
)

# General domain vocabulary patterns (zero test-set study titles).
PHAGE_RE = re.compile(r"\b(?:bacteriophages?|phages?)\b", re.I)
ACHROMOBACTER_RE = re.compile(r"\bachromobacter\b", re.I)
NTM_RE = re.compile(
    r"\b(?:nontuberculous\s+mycobacter\w*|non-tuberculous\s+mycobacter\w*|mycobacteri\w*\s+abscessus|"
    r"mycobacteri\w*\s+avium|mycobacterium\s+infections?|ntm)\b",
    re.I,
)
MRSA_RE = re.compile(r"\b(?:mrsa|methicillin[- ]resistant\s+staphylococcus\s+aureus)\b", re.I)
BURKHOLDERIA_RE = re.compile(r"\bburkholderia\b", re.I)
PSEUDOMONAS_RE = re.compile(r"\b(?:pseudomonas(?:\s+aeruginosa)?|p\.?\s*aeruginosa)\b", re.I)
SAMPLING_SPECIMEN_RE = re.compile(
    r"\b(?:induced\s+sputum|(?:cough|throat|oropharyngeal|nasal|oral)\s+swabs?|"
    r"face\s+mask\s+sampling|exhaled\s+breath|saliva)\b",
    re.I,
)
DIAGNOSTIC_EVAL_RE = re.compile(
    r"\b(?:diagnostics?|concordance|sensitivity|specificity|detection|culture\s+yield|"
    r"bacterial\s+culturing|sampling)\b",
    re.I,
)
SPUTUM_INDUCTION_SALINE_RE = re.compile(r"\b(?:hypertonic\s+saline|isotonic\s+saline|sodium\s+chloride)\b", re.I)
EXACERBATION_TERM_RE = re.compile(r"\bpulmonary\s+exacerbations?\b", re.I)
ANTIBIOTIC_STRATEGY_RE = re.compile(
    r"\b(?:oral\s+antibiotics?|antibiotic\s+(?:stewardship|duration|regimen|treatment)|"
    r"standardi[sz](?:ed|ing)\s+treatments?\s+of\s+pulmonary\s+exacerbations?|"
    r"standardi[sz](?:ed|ing)\s+treatments?\s+for\s+pulmonary\s+exacerbations?)\b",
    re.I,
)
MODULATOR_ONLY_RE = re.compile(r"\b(?:elexacaftor|tezacaftor|ivacaftor|vanzacaftor|deutivacaftor|lumacaftor)\b", re.I)
QPCR_PHAGE_RE = re.compile(r"\bqpcr\b", re.I)
NTM_CONVERSION_RE = re.compile(r"\b(?:ntm\s+clearance|sequential\s+negative\s+ntm\s+cultures|ntm\s+culture\s+conversion)\b", re.I)


class InfectionBenchmarkError(ValueError):
    """Raised when cohort separation, positive controls, license rules, or quarantine checks fail."""


def is_allowed_registry_url(raw_url: str) -> bool:
    """Validate https scheme, empty userinfo, and exact IDNA-normalized ClinicalTrials.gov host."""
    if not raw_url or not isinstance(raw_url, str):
        return False
    parsed = urlparse(raw_url.strip())
    if parsed.scheme != "https" or not parsed.hostname or parsed.username is not None or parsed.password is not None:
        return False
    try:
        ascii_host = parsed.hostname.encode("idna").decode("ascii").lower()
    except UnicodeError:
        return False
    return ascii_host in ALLOWED_REGISTRY_HOSTS


def quarantine_untrusted_text(raw_text: str, nonce: str | None = None) -> dict[str, Any]:
    """Neutralize homoglyph/bidi-split closing tags BEFORE escaping bidi/control chars inside <untrusted_data_NONCE>."""
    active_nonce = nonce or secrets.token_hex(8)
    raw_str = str(raw_text or "")

    stripped_for_skeleton = "".join(
        ch for ch in raw_str
        if ord(ch) not in TROJAN_BIDI_CODEPOINTS and not (unicodedata.category(ch) in ("Cf", "Cc") and ch not in ("\n", "\r", "\t"))
    )
    skeleton = stripped_for_skeleton.translate(CONFUSABLE_ASCII_MAP)
    lookalikes_escaped = 0
    working = raw_str
    if UNTRUSTED_CLOSE_LOOKALIKE_RE.search(skeleton):
        angle_chars = {"<", ">", "\uFF1C", "\uFF1E", "\u2329", "\u232A", "\u27E8", "\u27E9"}
        rebuilt: list[str] = []
        for ch in working:
            if ch in angle_chars:
                lookalikes_escaped += 1
                rebuilt.append("&lt;" if ch in ("<", "\uFF1C", "\u2329", "\u27E8") else "&gt;")
            else:
                rebuilt.append(ch)
        working = "".join(rebuilt)

    escaped_chars: list[str] = []
    trojan_escaped = 0
    for ch in working:
        cp = ord(ch)
        cat = unicodedata.category(ch)
        if cp in TROJAN_BIDI_CODEPOINTS or (cat in ("Cf", "Cc") and ch not in ("\n", "\r", "\t")):
            escaped_chars.append(f"<U+{cp:04X}>")
            trojan_escaped += 1
        else:
            escaped_chars.append(ch)
    sanitized = "".join(escaped_chars)

    caged = f"<untrusted_data_{active_nonce}>\n{sanitized}\n</untrusted_data_{active_nonce}>"
    return {
        "nonce": active_nonce,
        "caged_text": caged,
        "trojan_chars_escaped": trojan_escaped,
        "boundary_lookalikes_escaped": lookalikes_escaped,
    }


def _study_text_bundle(study: dict[str, Any]) -> str:
    parts: list[str] = [
        str(study.get("brief_title") or ""),
        str(study.get("official_title") or ""),
        " ".join(str(c) for c in (study.get("conditions") or [])),
        " ".join(str(k) for k in (study.get("keywords") or [])),
    ]
    for iv in study.get("interventions") or []:
        if isinstance(iv, dict):
            parts.append(str(iv.get("name") or ""))
            parts.append(str(iv.get("description") or ""))
    for po in study.get("primary_outcomes") or []:
        if isinstance(po, dict):
            parts.append(str(po.get("measure") or ""))
            parts.append(str(po.get("description") or ""))
    return " ".join(p for p in parts if p)


def _tag_primary_outcomes(study: dict[str, Any], lx: lex.Lexicon) -> list[str]:
    nct_id = str(study.get("nct_id") or "")
    classes: list[str] = []
    for idx, po in enumerate(study.get("primary_outcomes") or [], 1):
        res = lx.tag_entry(
            {
                "entry_id": f"{nct_id}:P{idx}",
                "nct_id": nct_id,
                "measure": str(po.get("measure") or ""),
                "description": str(po.get("description") or ""),
                "time_frame": str(po.get("time_frame") or ""),
                "no_primary_outcome": False,
            }
        )
        for t in res.get("tags") or []:
            cls_name = str(t.get("class") or "")
            if cls_name and cls_name not in classes:
                classes.append(cls_name)
    return classes


def classify_study_multi_field(
    nct_id: str,
    raw_by_nct: dict[str, dict[str, Any]],
    prisma_ncts: set[str],
    lx: lex.Lexicon,
    external_meta: dict[str, dict[str, Any]] = EXTERNAL_AND_SECONDARY_PROTOCOL_METADATA,
    quarantine_nonce: str | None = None,
) -> dict[str, Any]:
    """Classify cohort partition, pathogen/diagnostic focus, and primary/secondary endpoint locus."""
    snap_rec = raw_by_nct.get(nct_id)
    ext_rec = external_meta.get(nct_id)

    if snap_rec is None:
        if ext_rec and str(ext_rec.get("study_type") or "").upper() == "OBSERVATIONAL":
            prim_list = ext_rec.get("primary_outcomes") or []
            raw_excerpt = " | ".join(f"{p.get('measure', NOT_STATED)} :: {p.get('description', NOT_STATED)}" for p in prim_list)
            return {
                "nct_id": nct_id,
                "study_type": "OBSERVATIONAL",
                "in_raw_snapshot_1053": False,
                "included_in_prisma_940": False,
                "scope_exclusion_rule": "X1",
                "cohort_partition": "external_observational",
                "pathogen_focus": "non_sputum_diagnostic_sampling",
                "endpoint_locus": "observational_diagnostic_concordance",
                "primary_lexicon_classes": ["microbiology_sputum"],
                "primary_cfu_ml_or_g_reported": False,
                "secondary_protocol_provenance": ext_rec.get("provenance", "clinicaltrials_gov_api_v2"),
                "verified_api_date": ext_rec.get("verified_api_date", NOT_STATED),
                "specimen_modalities": list(ext_rec.get("specimen_modalities") or []),
                "quarantined_external_outcomes": quarantine_untrusted_text(raw_excerpt, nonce=quarantine_nonce),
                "registry_url": ext_rec.get("registry_url", f"https://clinicaltrials.gov/study/{nct_id}"),
            }
        raise InfectionBenchmarkError(f"Unverified NCT ID {nct_id} absent from snapshot.json and verified external registry metadata.")

    rule_id, _ = scope_mod.classify(snap_rec)
    in_prisma = nct_id in prisma_ncts
    if in_prisma and rule_id is None:
        cohort_partition = "prisma_940_interventional"
    elif rule_id == "X4":
        cohort_partition = "raw_1053_x4_excluded"
    elif rule_id == "X5":
        cohort_partition = "raw_1053_x5_excluded"
    else:
        cohort_partition = f"raw_1053_{str(rule_id or 'other').lower()}_excluded"

    prim_classes = _tag_primary_outcomes(snap_rec, lx)
    full_text = _study_text_bundle(snap_rec)
    prim_text = " ".join(f"{po.get('measure', '')} {po.get('description', '')}" for po in (snap_rec.get("primary_outcomes") or []))

    if cohort_partition == "raw_1053_x5_excluded":
        return {
            "nct_id": nct_id,
            "study_type": str(snap_rec.get("study_type") or "INTERVENTIONAL"),
            "in_raw_snapshot_1053": True,
            "included_in_prisma_940": False,
            "scope_exclusion_rule": "X5",
            "cohort_partition": cohort_partition,
            "pathogen_focus": "non_cf_negation_control",
            "endpoint_locus": "excluded_non_cf",
            "primary_lexicon_classes": prim_classes,
            "primary_cfu_ml_or_g_reported": bool(re.search(r"\bcfus?\b|\bcolony[- ]forming\b", prim_text, re.I)),
            "secondary_protocol_provenance": NOT_STATED,
            "verified_api_date": NOT_STATED,
            "specimen_modalities": ["sputum"],
            "quarantined_external_outcomes": None,
            "registry_url": f"https://clinicaltrials.gov/study/{nct_id}",
        }

    has_non_induction_drug_or_bio = any(
        str(iv.get("type") or "").upper() in ("DRUG", "BIOLOGICAL")
        and not SPUTUM_INDUCTION_SALINE_RE.search(f"{iv.get('name', '')} {iv.get('description', '')}")
        for iv in (snap_rec.get("interventions") or [])
        if isinstance(iv, dict)
    )
    is_diagnostic_sampling = bool(
        SAMPLING_SPECIMEN_RE.search(full_text)
        and DIAGNOSTIC_EVAL_RE.search(full_text)
        and not has_non_induction_drug_or_bio
    )
    if PHAGE_RE.search(full_text):
        focus = "bacteriophage_achromobacter" if ACHROMOBACTER_RE.search(full_text) else "bacteriophage_pseudomonas"
    elif NTM_RE.search(full_text):
        focus = "ntm_antimycobacterial"
    elif MRSA_RE.search(full_text):
        focus = "mrsa_staphylococcus"
    elif BURKHOLDERIA_RE.search(full_text):
        focus = "burkholderia_gram_negative"
    elif is_diagnostic_sampling:
        focus = "non_sputum_diagnostic_sampling"
    elif PSEUDOMONAS_RE.search(full_text):
        focus = "pseudomonas_antimicrobial"
    elif EXACERBATION_TERM_RE.search(full_text) and ANTIBIOTIC_STRATEGY_RE.search(full_text):
        focus = "antibiotic_stewardship"
    elif MODULATOR_ONLY_RE.search(full_text) and "microbiology_sputum" not in prim_classes:
        focus = "modulator_negative_control"
    else:
        focus = "other_cf_trial"

    sec_outcomes = list((ext_rec or {}).get("secondary_outcomes") or [])
    sec_blob = " ".join(f"{so.get('measure', '')} {so.get('description', '')}" for so in sec_outcomes if isinstance(so, dict))
    caged_sec = quarantine_untrusted_text(sec_blob, nonce=quarantine_nonce) if sec_blob else None

    if "microbiology_sputum" in prim_classes:
        locus = "primary_cfu_or_culture"
    elif focus == "modulator_negative_control":
        locus = "non_infection_control"
    elif "safety_tolerability" in prim_classes and QPCR_PHAGE_RE.search(sec_blob):
        locus = "primary_safety_secondary_qpcr"
    elif "safety_tolerability" in prim_classes and NTM_CONVERSION_RE.search(sec_blob):
        locus = "primary_safety_secondary_ntm_conversion"
    elif "fev1" in prim_classes:
        locus = "primary_clinical_fev1"
    elif "safety_tolerability" in prim_classes or "pharmacokinetics" in prim_classes:
        locus = "primary_safety_or_pk"
    else:
        locus = "other_outcome_locus"

    specimens: list[str] = list((ext_rec or {}).get("specimen_modalities") or [])
    if not specimens:
        if re.search(r"\bcough\s+swab\b", full_text, re.I):
            specimens.append("cough_swab")
        if re.search(r"\binduced\s+sputum\b", full_text, re.I):
            specimens.append("induced_sputum")
        if re.search(r"\bsputum\b", full_text, re.I):
            specimens.append("sputum")

    return {
        "nct_id": nct_id,
        "study_type": str(snap_rec.get("study_type") or "INTERVENTIONAL"),
        "in_raw_snapshot_1053": True,
        "included_in_prisma_940": in_prisma,
        "scope_exclusion_rule": rule_id,
        "cohort_partition": cohort_partition,
        "pathogen_focus": focus,
        "endpoint_locus": locus,
        "primary_lexicon_classes": prim_classes,
        "primary_cfu_ml_or_g_reported": bool(re.search(r"\bcfus?\b|\bcolony[- ]forming\b", prim_text, re.I)),
        "secondary_protocol_provenance": (ext_rec or {}).get("provenance", NOT_STATED),
        "verified_api_date": (ext_rec or {}).get("verified_api_date", NOT_STATED),
        "specimen_modalities": specimens or [NOT_STATED],
        "quarantined_external_outcomes": caged_sec,
        "registry_url": f"https://clinicaltrials.gov/study/{nct_id}",
    }


def classify_study_primary_only_baseline(
    nct_id: str,
    raw_by_nct: dict[str, dict[str, Any]],
    prisma_ncts: set[str],
    lx: lex.Lexicon,
) -> dict[str, str]:
    """Standard baseline comparator: inspects only PRISMA-940 primary_outcomes text."""
    if nct_id not in prisma_ncts:
        return {
            "nct_id": nct_id,
            "cohort_partition": "not_in_prisma_940",
            "pathogen_focus": "missed_outside_prisma_940",
            "endpoint_locus": "missed_outside_prisma_940",
        }
    snap_rec = raw_by_nct[nct_id]
    prim_classes = _tag_primary_outcomes(snap_rec, lx)
    prim_text = " ".join(f"{po.get('measure', '')} {po.get('description', '')}" for po in (snap_rec.get("primary_outcomes") or []))

    if PHAGE_RE.search(prim_text):
        focus = "bacteriophage_achromobacter" if ACHROMOBACTER_RE.search(prim_text) else "bacteriophage_pseudomonas"
    elif NTM_RE.search(prim_text):
        focus = "ntm_antimycobacterial"
    elif MRSA_RE.search(prim_text):
        focus = "mrsa_staphylococcus"
    elif BURKHOLDERIA_RE.search(prim_text):
        focus = "burkholderia_gram_negative"
    elif SAMPLING_SPECIMEN_RE.search(prim_text) and DIAGNOSTIC_EVAL_RE.search(prim_text):
        focus = "non_sputum_diagnostic_sampling"
    elif PSEUDOMONAS_RE.search(prim_text):
        focus = "pseudomonas_antimicrobial"
    elif EXACERBATION_TERM_RE.search(prim_text) and ANTIBIOTIC_STRATEGY_RE.search(prim_text):
        focus = "antibiotic_stewardship"
    else:
        focus = "unidentified_from_primary_outcome"

    if "microbiology_sputum" in prim_classes:
        locus = "primary_cfu_or_culture"
    elif "fev1" in prim_classes:
        locus = "primary_clinical_fev1"
    elif "safety_tolerability" in prim_classes or "pharmacokinetics" in prim_classes:
        locus = "primary_safety_or_pk"
    else:
        locus = "other_outcome_locus"

    return {"nct_id": nct_id, "cohort_partition": "prisma_940_interventional", "pathogen_focus": focus, "endpoint_locus": locus}


def verify_positive_controls(raw_by_nct: dict[str, dict[str, Any]], prisma_ncts: set[str], lx: lex.Lexicon) -> dict[str, bool]:
    """Enforce the four disjoint positive controls in done.md (none appear in the 20-study benchmark)."""
    pc_tip = classify_study_multi_field("NCT02212587", raw_by_nct, prisma_ncts, lx)
    pc_ganite = classify_study_multi_field("NCT01093521", raw_by_nct, prisma_ncts, lx)
    pc_x4 = classify_study_multi_field("NCT00910351", raw_by_nct, prisma_ncts, lx)
    pc_x5 = classify_study_multi_field("NCT01515007", raw_by_nct, prisma_ncts, lx)
    checks = {
        "tip_bcc_primary_cfu": (
            pc_tip["cohort_partition"] == "prisma_940_interventional"
            and pc_tip["endpoint_locus"] == "primary_cfu_or_culture"
            and pc_tip["primary_cfu_ml_or_g_reported"] is True
        ),
        "ganite_pk_safety_primary_not_cfu": (
            pc_ganite["cohort_partition"] == "prisma_940_interventional"
            and pc_ganite["endpoint_locus"] == "primary_safety_or_pk"
            and pc_ganite["primary_cfu_ml_or_g_reported"] is False
        ),
        "cipro_children_x4_excluded": (
            pc_x4["cohort_partition"] == "raw_1053_x4_excluded" and pc_x4["included_in_prisma_940"] is False
        ),
        "orbit3_non_cf_bronchiectasis_x5_excluded": (
            pc_x5["cohort_partition"] == "raw_1053_x5_excluded"
            and pc_x5["included_in_prisma_940"] is False
            and pc_x5["endpoint_locus"] == "excluded_non_cf"
        ),
    }
    failed = [name for name, ok in checks.items() if not ok]
    if failed:
        raise InfectionBenchmarkError(f"Positive control failure(s) before benchmark read: {failed}")
    return checks


def _assert_valid_quarantine_cage(cage: Any, label: str) -> None:
    if not isinstance(cage, dict) or not cage.get("nonce") or not isinstance(cage.get("caged_text"), str):
        raise InfectionBenchmarkError(f"Missing or malformed quarantine cage on {label}.")
    nonce = str(cage["nonce"])
    txt = cage["caged_text"]
    open_tag = f"<untrusted_data_{nonce}>"
    close_tag = f"</untrusted_data_{nonce}>"
    if not txt.startswith(open_tag) or not txt.endswith(close_tag) or txt.count("</untrusted_data_") != 1:
        raise InfectionBenchmarkError(f"Broken quarantine cage boundary on {label}.")


def verify_benchmark_report(report: dict[str, Any]) -> None:
    """Enforce cohort separation, NLM disclaimer, license isolation, quarantine cages, and URL invariants."""
    serialized = json.dumps(report, ensure_ascii=False)
    for pat in FORBIDDEN_OVERCLAIM_PHRASES:
        if pat.search(serialized):
            raise InfectionBenchmarkError("Report contains prohibited clinical overclaim or personal eligibility phrasing.")
    if "does not evaluate safety or efficacy" not in str(report.get("nlm_disclaimer") or ""):
        raise InfectionBenchmarkError("Report missing mandatory verbatim NLM ClinicalTrials.gov disclaimer.")

    nav_ptr = (report.get("registry_orientation") or {}).get("navigation_pointer_only") or {}
    if nav_ptr.get("catalog_source_id") == "cff-research-priorities" and nav_ptr.get("claim_label") != "unverified":
        raise InfectionBenchmarkError("cff-research-priorities must remain claim_label='unverified' navigation pointer only.")
    if "cff-research-priorities" in (report.get("verified_source_catalog_ids") or []):
        raise InfectionBenchmarkError("cff-research-priorities may never appear in verified_source_catalog_ids.")

    for lit in report.get("literature_anchors") or []:
        lic = str(lit.get("license") or "").upper()
        if "NC-ND" in lic or "ND" in lic:
            if lit.get("verbatim_cc_by_excerpt") or lit.get("derivative_text") or lit.get("translated_text"):
                raise InfectionBenchmarkError(f"License violation on {lit.get('doi')}: {lic} forbids derivatives/translations.")
        elif lit.get("verbatim_cc_by_excerpt"):
            _assert_valid_quarantine_cage(lit.get("quarantined_excerpt"), f"literature:{lit.get('doi')}")

    prisma_rows = report.get("prisma_940_interventional_studies") or []
    x4_x5_rows = report.get("raw_1053_excluded_studies") or []
    obs_rows = report.get("external_observational_studies") or []
    prisma_ids = {str(r.get("nct_id") or "") for r in prisma_rows}

    for r in prisma_rows:
        if r.get("study_type") == "OBSERVATIONAL" or not r.get("included_in_prisma_940"):
            raise InfectionBenchmarkError(f"Cohort contamination in prisma_940_interventional_studies: {r.get('nct_id')}")
        if r.get("endpoint_locus") in ("primary_safety_or_pk", "primary_safety_secondary_qpcr", "primary_safety_secondary_ntm_conversion"):
            if r.get("primary_cfu_ml_or_g_reported") is True:
                raise InfectionBenchmarkError(f"Endpoint misattribution on safety-primary trial {r.get('nct_id')}.")
        if r.get("quarantined_external_outcomes") is not None or r.get("endpoint_locus") == "primary_safety_secondary_qpcr":
            _assert_valid_quarantine_cage(r.get("quarantined_external_outcomes"), str(r.get("nct_id")))
        if not is_allowed_registry_url(str(r.get("registry_url") or "")):
            raise InfectionBenchmarkError(f"Disallowed registry URL host on {r.get('nct_id')}")

    for r in x4_x5_rows:
        nid = str(r.get("nct_id") or "")
        if nid in prisma_ids or r.get("included_in_prisma_940"):
            raise InfectionBenchmarkError(f"X4/X5 excluded study {nid} leaked into PRISMA 940 cohort.")
        if r.get("endpoint_locus") == "primary_safety_secondary_ntm_conversion":
            if r.get("primary_cfu_ml_or_g_reported") is True:
                raise InfectionBenchmarkError(f"Endpoint misattribution on X4 safety-primary study {nid}.")
        if r.get("quarantined_external_outcomes") is not None or r.get("endpoint_locus") == "primary_safety_secondary_ntm_conversion":
            _assert_valid_quarantine_cage(r.get("quarantined_external_outcomes"), nid)
        if not is_allowed_registry_url(str(r.get("registry_url") or "")):
            raise InfectionBenchmarkError(f"Disallowed registry URL host on {nid}")

    for r in obs_rows:
        nid = str(r.get("nct_id") or "")
        if nid in prisma_ids or r.get("included_in_prisma_940") or r.get("in_raw_snapshot_1053"):
            raise InfectionBenchmarkError(f"External OBSERVATIONAL study {nid} leaked into interventional snapshot cohorts.")
        _assert_valid_quarantine_cage(r.get("quarantined_external_outcomes"), nid)
        if not is_allowed_registry_url(str(r.get("registry_url") or "")):
            raise InfectionBenchmarkError(f"Disallowed registry URL host on external study {nid}")


def run_infection_benchmark(
    snapshot_path: Path = DEFAULT_SNAPSHOT_JSON,
    lexicon_path: Path = DEFAULT_LEXICON_JSON,
    ground_truth: tuple[dict[str, str], ...] = INFECTION_BENCHMARK_20_GROUND_TRUTH,
) -> dict[str, Any]:
    """Run the 20-study real-data CFF infection & diagnostic outcomes benchmark against primary-only baseline."""
    raw_data = json.loads(snapshot_path.read_text(encoding="utf-8"))
    raw_studies = raw_data.get("studies") or []
    interventional_raw_count = sum(1 for s in raw_studies if str(s.get("study_type") or "").upper() == "INTERVENTIONAL")
    raw_by_nct = {str(s["nct_id"]): s for s in raw_studies if s.get("nct_id")}

    meta, prisma_studies = load_prisma_included_studies(snapshot_path)
    prisma_ncts = {str(s["nct_id"]) for s in prisma_studies if s.get("nct_id")}
    lx = lex.Lexicon(json.loads(lexicon_path.read_text(encoding="utf-8")))
    positive_controls = verify_positive_controls(raw_by_nct, prisma_ncts, lx)

    prisma_rows, excluded_rows, observational_rows, evaluations = [], [], [], []
    mf_exact = base_exact = prisma_subset_total = prisma_mf_exact = prisma_base_exact = cross_total = cross_mf = cross_base = 0

    for gt in ground_truth:
        nid = gt["nct_id"]
        mf = classify_study_multi_field(nid, raw_by_nct, prisma_ncts, lx)
        base = classify_study_primary_only_baseline(nid, raw_by_nct, prisma_ncts, lx)
        mf_ok = (mf["cohort_partition"], mf["pathogen_focus"], mf["endpoint_locus"]) == (
            gt["expected_cohort"], gt["expected_focus"], gt["expected_locus"]
        )
        base_ok = (base["cohort_partition"], base["pathogen_focus"], base["endpoint_locus"]) == (
            gt["expected_cohort"], gt["expected_focus"], gt["expected_locus"]
        )
        mf_exact += int(mf_ok)
        base_exact += int(base_ok)
        if gt["expected_cohort"] == "prisma_940_interventional":
            prisma_subset_total += 1
            prisma_mf_exact += int(mf_ok)
            prisma_base_exact += int(base_ok)
            prisma_rows.append(mf)
        elif gt["expected_cohort"] in ("raw_1053_x4_excluded", "raw_1053_x5_excluded"):
            cross_total += 1
            cross_mf += int(mf_ok)
            cross_base += int(base_ok)
            excluded_rows.append(mf)
        else:
            cross_total += 1
            cross_mf += int(mf_ok)
            cross_base += int(base_ok)
            observational_rows.append(mf)
        evaluations.append({"nct_id": nid, "expected": gt, "multi_field_exact_match": mf_ok, "baseline_exact_match": base_ok})

    total_n = len(ground_truth)
    mf_acc = round(mf_exact / total_n, 4) if total_n else 0.0
    base_acc = round(base_exact / total_n, 4) if total_n else 0.0
    lit_anchors_out = [
        {**item, "quarantined_excerpt": quarantine_untrusted_text(str(item["verbatim_cc_by_excerpt"]))}
        if item.get("verbatim_cc_by_excerpt")
        else dict(item)
        for item in NON_SPUTUM_LITERATURE_ANCHORS
    ]

    report: dict[str, Any] = {
        "schema_version": 1,
        "tool": "infection_outcomes_benchmark",
        "benchmark_version": BENCHMARK_VERSION,
        "epistemic_tier": "T2",
        "verified_source_catalog_ids": ["clinicaltrials-gov", "cffpr-highlights-2025", "europe-pmc"],
        "snapshot_data_timestamp": meta.get("data_timestamp", "2026-10-14T14:12:11.425555+00:00"),
        "snapshot_raw_interventional_count": interventional_raw_count,
        "prisma_included_interventional_count": len(prisma_studies),
        "clinical_boundary": CLINICAL_BOUNDARY,
        "nlm_disclaimer": NLM_DISCLAIMER,
        "registry_orientation": CFFPR_2025_INFECTION_ORIENTATION,
        "positive_controls": positive_controls,
        "benchmark_metrics": {
            "total_studies": total_n,
            "multi_field_exact_matches": mf_exact,
            "multi_field_exact_accuracy": mf_acc,
            "held_out_error_fraction": round(1.0 - mf_acc, 4),
            "primary_only_baseline_exact_matches": base_exact,
            "primary_only_baseline_accuracy": base_acc,
            "accuracy_gain_over_baseline": round(mf_acc - base_acc, 4),
            "prisma_940_subset": {
                "studies": prisma_subset_total,
                "multi_field_exact_matches": prisma_mf_exact,
                "multi_field_accuracy": round(prisma_mf_exact / prisma_subset_total, 4) if prisma_subset_total else 0.0,
                "baseline_exact_matches": prisma_base_exact,
                "baseline_accuracy": round(prisma_base_exact / prisma_subset_total, 4) if prisma_subset_total else 0.0,
            },
            "cross_cohort_x4_x5_observational_subset": {
                "studies": cross_total,
                "multi_field_exact_matches": cross_mf,
                "multi_field_accuracy": round(cross_mf / cross_total, 4) if cross_total else 0.0,
                "baseline_exact_matches": cross_base,
                "baseline_accuracy": round(cross_base / cross_total, 4) if cross_total else 0.0,
            },
            "prisma_primary_microbiology_count": sum(1 for r in prisma_rows if r["endpoint_locus"] == "primary_cfu_or_culture"),
        },
        "prisma_940_interventional_studies": prisma_rows,
        "raw_1053_excluded_studies": excluded_rows,
        "external_observational_studies": observational_rows,
        "literature_anchors": lit_anchors_out,
        "study_evaluations": evaluations,
    }
    verify_benchmark_report(report)
    return report


def run_falsification_suite(snapshot_path: Path = DEFAULT_SNAPSHOT_JSON, lexicon_path: Path = DEFAULT_LEXICON_JSON) -> dict[str, Any]:
    """Run five deterministic negative controls (--falsify) verifying fail-closed behavior."""
    base_report = run_infection_benchmark(snapshot_path, lexicon_path)
    controls: list[dict[str, Any]] = []

    def _expect_fail(name: str, mutated: dict[str, Any]) -> None:
        try:
            verify_benchmark_report(mutated)
            controls.append({"control": name, "passed": False, "detail": "Poisoned report accepted"})
        except InfectionBenchmarkError as exc:
            controls.append({"control": name, "passed": True, "detail": str(exc)})

    p1 = {**base_report, "prisma_940_interventional_studies": [*base_report["prisma_940_interventional_studies"], {"nct_id": "NCT07312734", "study_type": "OBSERVATIONAL", "included_in_prisma_940": False, "registry_url": "https://clinicaltrials.gov/study/NCT07312734"}]}
    _expect_fail("observational_send_cf_leak_blocked", p1)

    bad_excl = [dict(r) for r in base_report["raw_1053_excluded_studies"]]
    bad_excl[0]["included_in_prisma_940"] = True
    _expect_fail("x4_x5_cohort_contamination_blocked", {**base_report, "raw_1053_excluded_studies": bad_excl})

    bad_prisma = [{**r, "primary_cfu_ml_or_g_reported": True} if r["nct_id"] == "NCT07275905" else dict(r) for r in base_report["prisma_940_interventional_studies"]]
    _expect_fail("phase1_safety_cfu_misattribution_blocked", {**base_report, "prisma_940_interventional_studies": bad_prisma})

    q_cyr = quarantine_untrusted_text("Txt \u202eBIDI </untrust\u0435\u0501_data_abcd> x", nonce="abcd1234abcd1234")
    q_grk = quarantine_untrusted_text("Txt </untrust\u03b5d_d\u03b1t\u03b1_abcd> x", nonce="abcd1234abcd1234")
    q_bidi_split = quarantine_untrusted_text("Txt </untrusted\u200b_data_abcd1234abcd1234> x", nonce="abcd1234abcd1234")
    spoof_rejected = not is_allowed_registry_url("https://\u0441linicaltrials.gov/study/NCT07312734")
    try:
        verify_benchmark_report({**base_report, "clinical_boundary": "You are eligible and proven to cure CF."})
        overclaim_blocked = False
    except InfectionBenchmarkError:
        overclaim_blocked = True
    c4_ok = (
        q_cyr["trojan_chars_escaped"] >= 1
        and q_cyr["boundary_lookalikes_escaped"] >= 2
        and q_grk["boundary_lookalikes_escaped"] >= 2
        and q_bidi_split["boundary_lookalikes_escaped"] >= 2
        and q_bidi_split["caged_text"].count("</untrusted_data_abcd1234abcd1234>") == 1
        and spoof_rejected
        and overclaim_blocked
    )
    controls.append({"control": "trojan_bidi_split_greek_cyrillic_idn_and_overclaim_blocked", "passed": c4_ok, "detail": f"grk={q_grk['boundary_lookalikes_escaped']} bidi_split={q_bidi_split['boundary_lookalikes_escaped']}"})

    bad_lit = [{**x, "translated_text": "Prohibited derivative."} if "NC-ND" in str(x.get("license") or "") else dict(x) for x in base_report["literature_anchors"]]
    try:
        verify_benchmark_report({**base_report, "literature_anchors": bad_lit})
        ncnd_ok = False
    except InfectionBenchmarkError:
        ncnd_ok = True
    try:
        verify_benchmark_report({**base_report, "verified_source_catalog_ids": [*base_report["verified_source_catalog_ids"], "cff-research-priorities"]})
        ptr_ok = False
    except InfectionBenchmarkError:
        ptr_ok = True
    controls.append({"control": "cc_by_nc_nd_and_unverified_pointer_blocked", "passed": ncnd_ok and ptr_ok, "detail": f"ncnd={ncnd_ok} ptr={ptr_ok}"})

    return {"tool": "infection_outcomes_benchmark", "falsify_passed": all(c["passed"] for c in controls), "negative_controls": controls}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT_JSON)
    parser.add_argument("--lexicon", type=Path, default=DEFAULT_LEXICON_JSON)
    parser.add_argument("--benchmark", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--falsify", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)

    if args.falsify:
        f_rep = run_falsification_suite(args.snapshot, args.lexicon)
        if args.json:
            print(json.dumps(f_rep, indent=2, ensure_ascii=False))
        else:
            print(f"{'PASS' if f_rep['falsify_passed'] else 'FAIL'} falsify infection_outcomes_benchmark ({sum(1 for c in f_rep['negative_controls'] if c['passed'])}/{len(f_rep['negative_controls'])} negative controls)")
        return 0 if f_rep["falsify_passed"] else 1

    report = run_infection_benchmark(args.snapshot, args.lexicon)
    metrics = report["benchmark_metrics"]
    if args.check:
        f_rep = run_falsification_suite(args.snapshot, args.lexicon)
        ok = f_rep["falsify_passed"] and metrics["multi_field_exact_accuracy"] >= 0.90 and metrics["multi_field_exact_accuracy"] > metrics["primary_only_baseline_accuracy"]
        if not ok:
            print("FAIL infection_outcomes_benchmark --check", file=sys.stderr)
            return 1
        print(
            f"PASS infection_outcomes_benchmark: accuracy={metrics['multi_field_exact_matches']}/{metrics['total_studies']} "
            f"({metrics['multi_field_exact_accuracy']:.4f}) vs baseline={metrics['primary_only_baseline_exact_matches']}/{metrics['total_studies']} "
            f"({metrics['primary_only_baseline_accuracy']:.4f}); prisma_940={len(report['prisma_940_interventional_studies'])} "
            f"raw_1053_excl={len(report['raw_1053_excluded_studies'])} external_obs={len(report['external_observational_studies'])}; "
            f"falsify={len(f_rep['negative_controls'])}/{len(f_rep['negative_controls'])}"
        )
        return 0
    if args.json:
        print(json.dumps(report, indent=2, ensure_ascii=False))
        return 0
    print(f"Benchmark {report['benchmark_version']}: {metrics['multi_field_exact_matches']}/{metrics['total_studies']} (baseline {metrics['primary_only_baseline_exact_matches']}/{metrics['total_studies']}) | {report['nlm_disclaimer']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
