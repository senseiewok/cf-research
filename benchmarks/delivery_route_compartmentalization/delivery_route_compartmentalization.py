#!/usr/bin/env python3
"""Delivery-Route Biomarker Compartmentalization Benchmark (`PRED-0015`).

Tests Engine B's pre-registered compartmentalization prediction (`benchmarks/delivery_route_compartmentalization/done.md`)
on measured public clinical-trial records in `tools/trial_atlas/snapshot.json` (`catalog.yaml` source
`clinicaltrials-gov`, `metadata.data_timestamp` `2026-11-23T02:01:59.090764+00:00` [parent Trial Atlas census `2026-10-09`], `940` PRISMA-included interventional CF trials):

- Arm 1 (`systemic_cftr`): Early-phase non-healthy-volunteer interventional CF trials testing
  systemic small-molecule CFTR modulators or systemic PTC readthrough agents (`interventions[]`).
- Arm 2 (`inhaled_cftr`): Early-phase non-healthy-volunteer interventional CF trials testing
  airway/respiratory-mucosal CFTR-restoring nucleic-acid, viral-vector, non-viral gene-transfer,
  splice-switching oligonucleotide, or topical-aminoglycoside PTC readthrough regimens (with
  investigational-route precedence over co-administered systemic modulators and a pure-route
  sensitivity check).

Invariants (`cf-research-context`, `discovery-practice`, `security-browsing`, `AH-1`..`AH-4`):
- 100% pure Python 3.10+ standard library + local `tools/trial_atlas/` modules; zero network/LLM calls.
- Strict `AH-1` verbatim NLM disclaimer and `AH-2` PRISMA-940 vs raw-1053 cohort boundary.
- Reads flattened `primary_outcomes[].measure` only (never assuming uncollected `secondary_outcomes`).
- Enforces CodeQL-safe identifiers (`RESEARCH_BOUNDARY`, `NLM_ATTRIBUTION`, `cl_transport`).
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import random
import re
import sys
from typing import Any

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
TRIAL_ATLAS_DIR = REPO_ROOT / "tools" / "trial_atlas"
if str(TRIAL_ATLAS_DIR) not in sys.path:
    sys.path.insert(0, str(TRIAL_ATLAS_DIR))

from consensus_audit import _to_raw_study  # noqa: E402
import scope as scope_mod  # noqa: E402
import snapshot as snap  # noqa: E402

BENCHMARK_VERSION = "2026-11-24.v1"
DEFAULT_SNAPSHOT_JSON = TRIAL_ATLAS_DIR / "snapshot.json"
DEFAULT_DONE_MD = HERE / "done.md"
DEFAULT_RESULTS_JSON = HERE / "results_2026-11-24.json"

NLM_ATTRIBUTION = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)
RESEARCH_BOUNDARY = (
    "Independent open-source research benchmark only (Tier T1/T2 registry metadata); "
    "not clinical guidance and not an assessment of personal trial eligibility or treatment efficacy."
)

EARLY_PHASES: frozenset[str] = frozenset({"EARLY_PHASE1", "PHASE1", "PHASE2"})

HEALTHY_VOLUNTEER_RE = re.compile(
    r"\bhealthy\s+(?:volunteer|volunteers|subject|subjects|participant|participants|adult|adults|male|males|female|females)\b",
    re.I,
)
CF_PATIENT_COHORT_RE = re.compile(
    r"\b(?:patients|subjects|participants|people|adults|children|adolescents|individuals)\s+with\s+cystic\s+fibrosis\b"
    r"|\bcystic\s+fibrosis\s+(?:patients|subjects|participants)\b",
    re.I,
)

# Named airway/respiratory-mucosal CFTR-restoring clinical candidates (`done.md` Section 2).
NAMED_INHALED_CFTR_RE = re.compile(
    r"\b(?:"
    r"vx[\s-]?522|mrt[\s-]?5005|arct[\s-]?032|rct[\s-]?2100|bmd[\s-]?003|lunar[\s-]?cf|"
    r"4d[\s-]?710|sp[\s-]?101|kb[\s-]?407|bi[\s-]?3720931|lenticlair|"
    r"pgm[\s-]?169|gl[\s-]?67a|"
    r"eluforsen|qr[\s-]?010|spl[\s-]?84(?:[\s-]?23[\s-]?1)?"
    r")\b",
    re.I,
)

CFTR_GENETIC_TARGET_RE = re.compile(
    r"\b(?:"
    r"cftr\s+(?:mrna|cdna|dna|plasmid|gene|vector|antisense|oligonucleotide)|"
    r"(?:mrna|cdna|plasmid|gene\s+transfer|gene\s+therapy|viral\s+vector|adenovir\w*|ad5|adgv\w*|raav\w*|aav\w*|lentivir\w*|liposom\w*)"
    r"(?:[^.;]{0,60})\bcftr\b|"
    r"\bcftr\b(?:[^.;]{0,60})(?:mrna|cdna|plasmid|gene\s+transfer|gene\s+therapy|adenovir\w*|ad5|adgv\w*|raav\w*|aav\w*|lentivir\w*|liposom\w*|antisense\s+oligonucleotide)|"
    r"ad[\s-]?cftr|aav2?[\s-]?cftr|"
    r"antisense\s+oligonucleotide(?:[^.;]{0,50})(?:cftr|f508del|3849\+10kb)"
    r")",
    re.I,
)

RESPIRATORY_MUCOSAL_ROUTE_RE = re.compile(
    r"\b(?:"
    r"inhal\w*|nebuliz\w*|nebulis\w*|aerosol\w*|pulmonary|endobronchial|bronchoscop\w*|"
    r"intranasal|nasal|maxillary\s+sinus|sinus|airway"
    r")\b",
    re.I,
)

TOPICAL_AMINOGLYCOSIDE_RE = re.compile(r"\b(?:gentamicin|amikacin)\b", re.I)
PTC_READTHROUGH_INTENT_RE = re.compile(
    r"\b(?:nonsense|stop\s+codon|premature\s+stop|premature\s+termination|readthrough|read-through|w1282x|g542x)\b",
    re.I,
)

# Canonical systemic CFTR modulators and systemic PTC readthrough agents (`done.md` Section 2).
SYSTEMIC_CFTR_AGENT_RE = re.compile(
    r"\b(?:"
    r"ivacaftor|iva|vx[\s-]?770|kalydeco|"
    r"deutivacaftor|vx[\s-]?561|ctp[\s-]?656|"
    r"lumacaftor|lum/iva|vx[\s-]?809|orkambi|"
    r"tezacaftor|tez/iva|vx[\s-]?661|symdeko|symkevi|"
    r"elexacaftor|elx/tez/iva|vx[\s-]?445|trikafta|kaftrio|"
    r"vanzacaftor|vx[\s-]?121|alyftrek|"
    r"vx[\s-]?659|vx[\s-]?440|vx[\s-]?152|vx[\s-]?581|"
    r"galicaftor|abbv[\s-]?2222|glpg[\s-]?2222|"
    r"navocaftor|abbv[\s-]?3067|glpg[\s-]?3067|"
    r"glpg[\s-]?1837|glpg[\s-]?2451|glpg[\s-]?2737|"
    r"abbv[\s-]?191|abbv[\s-]?576|"
    r"posenacaftor|pti[\s-]?801|nesolicaftor|pti[\s-]?428|dirocaftor|pti[\s-]?808|"
    r"fdl[\s-]?169|icenticaftor|qbw[\s-]?251|cavosonstat|n[\s-]?9115|"
    r"4[\s-]?phenylbutyrate|sodium\s+phenylbutyrate|buphenyl|"
    r"ataluren|ptc[\s-]?124|elx[\s-]?02"
    r")\b",
    re.I,
)

# Systemic CFTR modulators specifically (for dual-route co-administration flag on inhaled_cftr).
SYSTEMIC_MODULATOR_ONLY_RE = re.compile(
    r"\b(?:"
    r"ivacaftor|iva|vx[\s-]?770|kalydeco|deutivacaftor|vx[\s-]?561|ctp[\s-]?656|"
    r"lumacaftor|lum/iva|vx[\s-]?809|orkambi|tezacaftor|tez/iva|vx[\s-]?661|symdeko|symkevi|"
    r"elexacaftor|elx/tez/iva|vx[\s-]?445|trikafta|kaftrio|vanzacaftor|vx[\s-]?121|alyftrek|"
    r"vx[\s-]?659|vx[\s-]?440|vx[\s-]?152|vx[\s-]?581|"
    r"galicaftor|abbv[\s-]?2222|glpg[\s-]?2222|navocaftor|abbv[\s-]?3067|glpg[\s-]?3067|"
    r"glpg[\s-]?1837|glpg[\s-]?2451|glpg[\s-]?2737|abbv[\s-]?191|abbv[\s-]?576|"
    r"posenacaftor|pti[\s-]?801|nesolicaftor|pti[\s-]?428|dirocaftor|pti[\s-]?808|"
    r"fdl[\s-]?169|icenticaftor|qbw[\s-]?251|cavosonstat|n[\s-]?9115|"
    r"4[\s-]?phenylbutyrate|sodium\s+phenylbutyrate|buphenyl"
    r")\b",
    re.I,
)

WITHDRAWAL_OR_DDI_TITLE_RE = re.compile(
    r"\b(?:"
    r"discontinuation|stopping|withdrawal\s+of\s+hypertonic|simplify|"
    r"drug[\s-]drug\s+interaction|midazolam|oral\s+contraceptive|itraconazole|rifampin"
    r")\b",
    re.I,
)

NON_CFTR_TOPICAL_AGENT_RE = re.compile(
    r"\b(?:"
    r"etd[\s-]?001|aro[\s-]?enac|vx[\s-]?371|p[\s-]?1037|spx[\s-]?101|bi[\s-]?1265162|"
    r"azd[\s-]?5634|qbw[\s-]?276|amiloride|denufosol|"
    r"hypertonic\s+saline|dornase\s+alfa|mannitol|"
    r"bx[\s-]?004|ap[\s-]?pa02|achromophage|tobramycin|aztreonam|colistin|gallium"
    r")\b",
    re.I,
)

SWEAT_CHLORIDE_MEASURE_RE = re.compile(
    r"\bsweat\s*(?:\[\s*)?(?:chloride|cl\b|cl-)"
    r"|\bchloride\s+concentration\s+in\s+sweat\b",
    re.I,
)
DIAGNOSTIC_EXCLUSION_RE = re.compile(
    r"\b(?:diagnosis|diagnostic|inclusion|eligibility|screening\s+criteria)\b",
    re.I,
)

AIRWAY_OR_SAFETY_MEASURE_RE = re.compile(
    r"\b(?:"
    r"fev1|fev\s*1|ppfev1|forced\s+expiratory\s+volume|spirometr\w*|"
    r"lung\s+clearance\s+index|lci(?:\s*2\.5)?|multiple\s+breath\s+washout|"
    r"nasal\s+potential\s+difference|npd|potential\s+difference|mucociliary|"
    r"adverse\s+event\w*|teae\w*|sae\w*|safety|tolerability|dose[\s-]limiting|toxicit\w*|vital\s+signs|laboratory"
    r")\b",
    re.I,
)

ALLOWED_INTERVENTION_TYPES: frozenset[str] = frozenset({"DRUG", "BIOLOGICAL", "GENETIC"})
REQUIRED_OUTCOME_KEYS: frozenset[str] = frozenset({"measure", "description", "time_frame"})
Z_95 = 1.959963984540054


class CompartmentalizationBenchmarkError(RuntimeError):
    """Raised when schema, positive-control, or falsification invariants fail."""


def assert_flattened_study_schema(study: dict[str, Any]) -> None:
    """Enforce `done.md` Section 1 flattened schema check on every PRISMA-included study."""
    nct_id = study.get("nct_id")
    if "outcomesModule" in study or "secondary_outcomes" in study:
        raise CompartmentalizationBenchmarkError(
            f"Study {nct_id} carries unflattened outcomesModule or uncollected secondary_outcomes key"
        )
    if not isinstance(study.get("phases"), list):
        raise CompartmentalizationBenchmarkError(f"Study {nct_id} missing flattened phases list")
    if not isinstance(study.get("interventions"), list):
        raise CompartmentalizationBenchmarkError(f"Study {nct_id} missing flattened interventions list")
    if not isinstance(study.get("brief_title"), str) or not isinstance(study.get("official_title"), str):
        raise CompartmentalizationBenchmarkError(f"Study {nct_id} missing flattened brief_title/official_title strings")
    outcomes = study.get("primary_outcomes")
    if not isinstance(outcomes, list):
        raise CompartmentalizationBenchmarkError(f"Study {nct_id} missing flattened primary_outcomes list")
    for idx, item in enumerate(outcomes):
        if not isinstance(item, dict) or not REQUIRED_OUTCOME_KEYS.issubset(item.keys()):
            raise CompartmentalizationBenchmarkError(
                f"Study {nct_id} primary_outcomes[{idx}] lacks required flattened keys {sorted(REQUIRED_OUTCOME_KEYS)}"
            )


def sha256_lf_text(path: Path) -> str:
    text = path.read_text(encoding="utf-8").replace("\r\n", "\n")
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def is_early_phase(phases: list[str] | None) -> bool:
    if not phases:
        return False
    norm = {str(p).strip().upper() for p in phases if str(p).strip()}
    return bool(norm) and norm <= EARLY_PHASES


def is_healthy_volunteer_only(brief_title: str, official_title: str) -> bool:
    title_corpus = f"{brief_title or ''} {official_title or ''}"
    if not HEALTHY_VOLUNTEER_RE.search(title_corpus):
        return False
    return not bool(CF_PATIENT_COHORT_RE.search(title_corpus))


def _extract_interventions(study: dict[str, Any]) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for raw in study.get("interventions") or []:
        if isinstance(raw, dict):
            out.append((str(raw.get("name") or ""), str(raw.get("type") or "").strip().upper()))
        elif isinstance(raw, str):
            out.append((raw, ""))
    return out


def classify_study_regimen(study: dict[str, Any]) -> dict[str, Any]:
    """Deterministically classify a study into `inhaled_cftr`, `systemic_cftr`, or `non_cftr_or_other`."""
    brief_title = str(study.get("brief_title") or "")
    official_title = str(study.get("official_title") or "")
    keywords = " ".join(str(k) for k in (study.get("keywords") or []))
    interventions = _extract_interventions(study)
    int_names_all = " ".join(name for name, _ in interventions)
    active_int_names = " ".join(
        name for name, itype in interventions if itype in ALLOWED_INTERVENTION_TYPES
    )
    regimen_corpus = f"{brief_title} {official_title} {int_names_all} {keywords}"
    title_corpus = f"{brief_title} {official_title}"

    has_sys_mod_intervention = bool(SYSTEMIC_MODULATOR_ONLY_RE.search(active_int_names))
    has_non_cftr_topical = bool(NON_CFTR_TOPICAL_AGENT_RE.search(f"{title_corpus} {active_int_names}"))

    # 1. Investigational-route precedence: Arm 2 (`inhaled_cftr`) evaluated first.
    named_inhaled = bool(NAMED_INHALED_CFTR_RE.search(regimen_corpus))
    genetic_respiratory = bool(
        not has_non_cftr_topical
        and CFTR_GENETIC_TARGET_RE.search(regimen_corpus)
        and RESPIRATORY_MUCOSAL_ROUTE_RE.search(regimen_corpus)
    )
    topical_ptc_aminoglycoside = bool(
        not has_non_cftr_topical
        and TOPICAL_AMINOGLYCOSIDE_RE.search(regimen_corpus)
        and RESPIRATORY_MUCOSAL_ROUTE_RE.search(regimen_corpus)
        and PTC_READTHROUGH_INTENT_RE.search(regimen_corpus)
    )
    if named_inhaled or genetic_respiratory or topical_ptc_aminoglycoside:
        return {
            "arm": "inhaled_cftr",
            "dual_route_modulator_co_admin": has_sys_mod_intervention,
            "rule_detail": (
                "named_inhaled_cftr"
                if named_inhaled
                else ("genetic_respiratory_cftr" if genetic_respiratory else "topical_ptc_aminoglycoside")
            ),
        }

    # 2. Strict non-CFTR topical exclusion (`done.md` Section 2) before Arm 1 (`systemic_cftr`).
    if has_non_cftr_topical:
        return {
            "arm": "non_cftr_or_other",
            "dual_route_modulator_co_admin": False,
            "rule_detail": "excluded_non_cftr_topical",
        }

    # 3. Arm 1 (`systemic_cftr`): active DRUG/BIOLOGICAL/GENETIC intervention matches canonical systemic CFTR list.
    if SYSTEMIC_CFTR_AGENT_RE.search(active_int_names):
        if WITHDRAWAL_OR_DDI_TITLE_RE.search(title_corpus):
            return {
                "arm": "non_cftr_or_other",
                "dual_route_modulator_co_admin": False,
                "rule_detail": "excluded_withdrawal_or_ddi",
            }
        return {
            "arm": "systemic_cftr",
            "dual_route_modulator_co_admin": False,
            "rule_detail": "systemic_cftr_intervention",
        }

    return {
        "arm": "non_cftr_or_other",
        "dual_route_modulator_co_admin": False,
        "rule_detail": "other_regimen",
    }


def classify_primary_outcomes(study: dict[str, Any]) -> dict[str, Any]:
    outcomes = study.get("primary_outcomes")
    if not isinstance(outcomes, list):
        raise CompartmentalizationBenchmarkError(
            f"Study {study.get('nct_id')} missing flattened primary_outcomes list"
        )
    has_sweat_cl = False
    has_airway_or_safety = False
    matched_sweat_measures: list[str] = []
    for idx, o in enumerate(outcomes):
        if not isinstance(o, dict) or not REQUIRED_OUTCOME_KEYS.issubset(o.keys()):
            raise CompartmentalizationBenchmarkError(
                f"Study {study.get('nct_id')} primary_outcomes[{idx}] missing flattened keys"
            )
        measure = str(o.get("measure") or "").strip()
        if not measure:
            continue
        if SWEAT_CHLORIDE_MEASURE_RE.search(measure) and not DIAGNOSTIC_EXCLUSION_RE.search(measure):
            has_sweat_cl = True
            matched_sweat_measures.append(measure)
        if AIRWAY_OR_SAFETY_MEASURE_RE.search(measure):
            has_airway_or_safety = True
    return {
        "has_primary_sweat_chloride": has_sweat_cl,
        "has_primary_airway_or_safety": has_airway_or_safety,
        "matched_sweat_measures": matched_sweat_measures,
    }


def fisher_hypergeometric_exact(n1: int, k1: int, n2: int, k2: int) -> dict[str, Any]:
    """Compute exact hypergeometric one-sided (upper tail for k1) and two-sided p-values in Q."""
    if n1 < 0 or n2 < 0 or not (0 <= k1 <= n1) or not (0 <= k2 <= n2):
        raise ValueError(f"Invalid 2x2 counts: n1={n1}, k1={k1}, n2={n2}, k2={k2}")
    total_n = n1 + n2
    total_k = k1 + k2
    denom = math.comb(total_n, total_k)
    j_min = max(0, total_k - n2)
    j_max = min(n1, total_k)

    pmf: dict[int, Fraction] = {}
    for j in range(j_min, j_max + 1):
        pmf[j] = Fraction(math.comb(n1, j) * math.comb(n2, total_k - j), denom)

    support_sum = sum(pmf.values(), Fraction(0, 1))
    if support_sum != Fraction(1, 1):
        raise CompartmentalizationBenchmarkError(f"Hypergeometric support sum != 1: {support_sum}")

    p_one_sided = sum((pmf[j] for j in range(k1, j_max + 1)), Fraction(0, 1))
    p_obs = pmf[k1]
    p_two_sided = sum((prob for prob in pmf.values() if prob <= p_obs), Fraction(0, 1))
    return {
        "p_one_sided_exact": f"{p_one_sided.numerator}/{p_one_sided.denominator}",
        "p_one_sided_float": float(p_one_sided),
        "p_one_sided_fraction": p_one_sided,
        "p_two_sided_exact": f"{p_two_sided.numerator}/{p_two_sided.denominator}",
        "p_two_sided_float": float(p_two_sided),
        "support_sum_exact": "1/1",
    }


def _wilson_single_ci(k: int, n: int, z: float = Z_95) -> tuple[float, float]:
    if n <= 0:
        raise ValueError("Sample size n must be positive")
    z2 = z * z
    p_hat = k / n
    denom = 2.0 * (n + z2)
    disc = math.sqrt(z2 + 4.0 * k * (1.0 - p_hat))
    low = 0.0 if k == 0 else (2.0 * k + z2 - z * disc) / denom
    high = 1.0 if k == n else (2.0 * k + z2 + z * disc) / denom
    return low, high


def newcombe_wilson_rd_ci95(n1: int, k1: int, n2: int, k2: int, z: float = Z_95) -> tuple[float, float]:
    p1 = k1 / n1
    p2 = k2 / n2
    rd = p1 - p2
    l1, u1 = _wilson_single_ci(k1, n1, z)
    l2, u2 = _wilson_single_ci(k2, n2, z)
    low = rd - math.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2)
    high = rd + math.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)
    return max(-1.0, low), min(1.0, high)


def _mn_restricted_p2(x1: int, n1: int, x2: int, n2: int, delta: float) -> float:
    """Miettinen-Nurminen (1985) restricted MLE of p2 under p1 - p2 = delta."""
    low_b = max(0.0, -delta)
    high_b = min(1.0, 1.0 - delta)
    if high_b <= low_b + 1e-15:
        return 0.5 * (low_b + high_b)

    total_n = float(n1 + n2)
    total_x = float(x1 + x2)
    l3 = total_n
    l2 = (n1 + 2 * n2) * delta - total_n - total_x
    l1 = (n2 * delta - total_n - 2 * x2) * delta + total_x
    l0 = x2 * delta * (1.0 - delta)

    q = (l2 * l2) / (9.0 * l3 * l3) - l1 / (3.0 * l3)
    r = (l2 ** 3) / (27.0 * (l3 ** 3)) - (l1 * l2) / (6.0 * l3 * l3) + l0 / (2.0 * l3)
    if q <= 1e-30:
        cand = -l2 / (3.0 * l3)
    else:
        u = math.copysign(math.sqrt(q), r if r != 0.0 else 1.0)
        cos_arg = max(-1.0, min(1.0, r / (u ** 3)))
        w = (math.pi + math.acos(cos_arg)) / 3.0
        cand = 2.0 * u * math.cos(w) - l2 / (3.0 * l3)
    return max(low_b, min(high_b, cand))


def _mn_score_z(x1: int, n1: int, x2: int, n2: int, delta: float) -> float:
    num = (x1 / n1) - (x2 / n2) - delta
    if abs(num) < 1e-15:
        return 0.0
    p2_tilde = _mn_restricted_p2(x1, n1, x2, n2, delta)
    p1_tilde = min(1.0, max(0.0, p2_tilde + delta))
    total_n = n1 + n2
    var = (p1_tilde * (1.0 - p1_tilde) / n1 + p2_tilde * (1.0 - p2_tilde) / n2) * (total_n / (total_n - 1.0))
    if var <= 1e-24:
        return math.copysign(1e12, num)
    return num / math.sqrt(var)


def chan_zhang_unconditional_rd_ci95(
    n1: int, k1: int, n2: int, k2: int, alpha_two_sided: float = 0.05
) -> tuple[float, float]:
    """Chan & Zhang (1999) exact unconditional profile-score two-sided CI for p1 - p2."""
    alpha_tail = 0.5 * alpha_two_sided
    comb1 = [math.comb(n1, i) for i in range(n1 + 1)]
    comb2 = [math.comb(n2, j) for j in range(n2 + 1)]

    def binom_pmf(n: int, combs: list[int], p: float) -> list[float]:
        if p <= 0.0:
            return [1.0] + [0.0] * n
        if p >= 1.0:
            return [0.0] * n + [1.0]
        q = 1.0 - p
        return [combs[i] * (p ** i) * (q ** (n - i)) for i in range(n + 1)]

    def tail_supremum(delta: float, upper: bool) -> float:
        low_p2 = max(0.0, -delta)
        high_p2 = min(1.0, 1.0 - delta)
        z_obs = _mn_score_z(k1, n1, k2, n2, delta)

        active_pairs: list[tuple[int, int]] = []
        for x1 in range(n1 + 1):
            for x2 in range(n2 + 1):
                z_xy = _mn_score_z(x1, n1, x2, n2, delta)
                if upper:
                    if z_xy >= z_obs - 1e-11:
                        active_pairs.append((x1, x2))
                else:
                    if z_xy <= z_obs + 1e-11:
                        active_pairs.append((x1, x2))

        if not active_pairs:
            return 0.0
        if len(active_pairs) == (n1 + 1) * (n2 + 1):
            return 1.0

        def prob_at(p2: float) -> float:
            p1 = min(1.0, max(0.0, p2 + delta))
            b1 = binom_pmf(n1, comb1, p1)
            b2 = binom_pmf(n2, comb2, p2)
            total = 0.0
            for x1, x2 in active_pairs:
                total += b1[x1] * b2[x2]
            return total

        if high_p2 <= low_p2 + 1e-14:
            return prob_at(0.5 * (low_p2 + high_p2))

        grid_n = 120
        best_p2 = low_p2
        best_val = -1.0
        for idx in range(grid_n + 1):
            p2 = low_p2 + (high_p2 - low_p2) * (idx / grid_n)
            val = prob_at(p2)
            if val > best_val:
                best_val = val
                best_p2 = p2

        # Local golden-section refinement around best grid point.
        step = (high_p2 - low_p2) / grid_n
        a = max(low_p2, best_p2 - step)
        b = min(high_p2, best_p2 + step)
        inv_phi = (math.sqrt(5.0) - 1.0) / 2.0
        for _ in range(24):
            c = b - inv_phi * (b - a)
            d = a + inv_phi * (b - a)
            fc = prob_at(c)
            fd = prob_at(d)
            if fc > best_val:
                best_val = fc
            if fd > best_val:
                best_val = fd
            if fc < fd:
                a = c
            else:
                b = d
        return min(1.0, max(0.0, best_val))

    rd_obs = (k1 / n1) - (k2 / n2)

    # Lower bound: smallest delta in [-1, rd_obs] with upper-tail supremum > alpha_tail.
    if k1 == 0 and k2 == n2:
        low_ci = -1.0
    else:
        lo, hi = -1.0 + 1e-12, rd_obs
        if tail_supremum(lo, upper=True) > alpha_tail:
            low_ci = -1.0
        else:
            for _ in range(42):
                mid = 0.5 * (lo + hi)
                if tail_supremum(mid, upper=True) > alpha_tail:
                    hi = mid
                else:
                    lo = mid
            low_ci = hi

    # Upper bound: largest delta in [rd_obs, 1] with lower-tail supremum > alpha_tail.
    if k1 == n1 and k2 == 0:
        high_ci = 1.0
    else:
        lo, hi = rd_obs, 1.0 - 1e-12
        if tail_supremum(hi, upper=False) > alpha_tail:
            high_ci = 1.0
        else:
            for _ in range(42):
                mid = 0.5 * (lo + hi)
                if tail_supremum(mid, upper=False) > alpha_tail:
                    lo = mid
                else:
                    hi = mid
            high_ci = lo

    return low_ci, high_ci


def verify_positive_controls(prisma_by_nct: dict[str, dict[str, Any]], early_inhaled_rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Run all three pre-registered positive controls before reporting observed proportions."""
    # Control 1: Synthetic fixture (12, 9) vs (10, 0)
    cal_fisher = fisher_hypergeometric_exact(12, 9, 10, 0)
    if cal_fisher["p_one_sided_fraction"] != Fraction(11, 24871):
        raise CompartmentalizationBenchmarkError(
            f"Positive control 1 failed: expected 11/24871, got {cal_fisher['p_one_sided_exact']}"
        )
    nw_lo, _ = newcombe_wilson_rd_ci95(12, 9, 10, 0)
    cz_lo, _ = chan_zhang_unconditional_rd_ci95(12, 9, 10, 0)
    if not (nw_lo > 0.35 and cz_lo > 0.35):
        raise CompartmentalizationBenchmarkError(
            f"Positive control 1 CI failed: nw_lo={nw_lo:.6f}, cz_lo={cz_lo:.6f}"
        )

    # Control 2: Known trial regimen-route controls
    expected_controls = (
        ("NCT05668741", "inhaled_cftr", True),
        ("NCT05248230", "inhaled_cftr", False),
        ("NCT06429176", "inhaled_cftr", False),
        ("NCT00846781", "non_cftr_or_other", False),
    )
    verified_ids: dict[str, dict[str, Any]] = {}
    for nct_id, exp_arm, exp_dual in expected_controls:
        study = prisma_by_nct.get(nct_id)
        if study is None:
            raise CompartmentalizationBenchmarkError(f"Positive control 2 missing {nct_id} in PRISMA-940")
        cls_res = classify_study_regimen(study)
        if cls_res["arm"] != exp_arm or bool(cls_res["dual_route_modulator_co_admin"]) != exp_dual:
            raise CompartmentalizationBenchmarkError(
                f"Positive control 2 mismatch on {nct_id}: got {cls_res}, expected arm={exp_arm}, dual={exp_dual}"
            )
        verified_ids[nct_id] = cls_res

    # Control 3: Airway/safety primary-endpoint specificity on early-phase inhaled_cftr (>= 80%)
    if not early_inhaled_rows:
        raise CompartmentalizationBenchmarkError("Positive control 3 failed: zero early-phase inhaled_cftr studies")
    airway_or_safety_n = sum(1 for r in early_inhaled_rows if r["has_primary_airway_or_safety"])
    specificity_ratio = airway_or_safety_n / len(early_inhaled_rows)
    if specificity_ratio < 0.80:
        raise CompartmentalizationBenchmarkError(
            f"Positive control 3 failed: inhaled airway/safety specificity {airway_or_safety_n}/{len(early_inhaled_rows)} = {specificity_ratio:.4f} < 0.80"
        )

    return {
        "calibration_fixture_12_9_10_0": {
            "p_fisher_1s_exact": cal_fisher["p_one_sided_exact"],
            "newcombe_wilson_ci95_low": round(nw_lo, 6),
            "chan_zhang_ci95_low": round(cz_lo, 6),
        },
        "known_trial_controls": verified_ids,
        "inhaled_airway_or_safety_specificity": {
            "matched": airway_or_safety_n,
            "total": len(early_inhaled_rows),
            "proportion": round(specificity_ratio, 6),
        },
    }


def evaluate_snapshot(
    snapshot_path: Path = DEFAULT_SNAPSHOT_JSON,
    done_path: Path = DEFAULT_DONE_MD,
) -> dict[str, Any]:
    raw_doc = json.loads(snapshot_path.read_text(encoding="utf-8"))
    meta = raw_doc.get("metadata") or {}
    studies = raw_doc.get("studies") or []
    raw_interventional_n = len(studies)

    prisma_studies: list[dict[str, Any]] = []
    excluded_counts = {"X1": 0, "X2": 0, "X3": 0, "X4": 0, "X5": 0}
    for s in studies:
        rec = snap.study_record(_to_raw_study(s))
        rule, _ = scope_mod.classify(rec)
        if rule is None:
            assert_flattened_study_schema(s)
            prisma_studies.append(s)
        elif rule in excluded_counts:
            excluded_counts[rule] += 1

    if raw_interventional_n != 1053 or len(prisma_studies) != 940:
        raise CompartmentalizationBenchmarkError(
            f"AH-2 cohort mismatch: raw={raw_interventional_n} (expected 1053), prisma={len(prisma_studies)} (expected 940)"
        )

    prisma_by_nct = {str(s["nct_id"]): s for s in prisma_studies}

    early_systemic: list[dict[str, Any]] = []
    early_inhaled: list[dict[str, Any]] = []
    early_hv_excluded_n = 0
    all_phase_systemic: list[dict[str, Any]] = []
    all_phase_inhaled: list[dict[str, Any]] = []

    for s in prisma_studies:
        reg = classify_study_regimen(s)
        out = classify_primary_outcomes(s)
        hv_only = is_healthy_volunteer_only(str(s.get("brief_title") or ""), str(s.get("official_title") or ""))
        early = is_early_phase(s.get("phases"))

        row = {
            "nct_id": str(s["nct_id"]),
            "brief_title": str(s.get("brief_title") or ""),
            "phases": list(s.get("phases") or []),
            "overall_status": str(s.get("overall_status") or ""),
            "arm": reg["arm"],
            "dual_route_modulator_co_admin": reg["dual_route_modulator_co_admin"],
            "rule_detail": reg["rule_detail"],
            "has_primary_sweat_chloride": out["has_primary_sweat_chloride"],
            "has_primary_airway_or_safety": out["has_primary_airway_or_safety"],
            "matched_sweat_measures": out["matched_sweat_measures"],
        }

        if not hv_only:
            if reg["arm"] == "systemic_cftr":
                all_phase_systemic.append(row)
            elif reg["arm"] == "inhaled_cftr":
                all_phase_inhaled.append(row)

        if early:
            if hv_only:
                early_hv_excluded_n += 1
                continue
            if reg["arm"] == "systemic_cftr":
                early_systemic.append(row)
            elif reg["arm"] == "inhaled_cftr":
                early_inhaled.append(row)

    # Verify positive controls BEFORE computing outcome proportions!
    pos_controls = verify_positive_controls(prisma_by_nct, early_inhaled)

    n_sys = len(early_systemic)
    k_sys = sum(1 for r in early_systemic if r["has_primary_sweat_chloride"])
    n_inh = len(early_inhaled)
    k_inh = sum(1 for r in early_inhaled if r["has_primary_sweat_chloride"])

    pure_inhaled = [r for r in early_inhaled if not r["dual_route_modulator_co_admin"]]
    n_inh_pure = len(pure_inhaled)
    k_inh_pure = sum(1 for r in pure_inhaled if r["has_primary_sweat_chloride"])

    p_sys = k_sys / n_sys if n_sys else 0.0
    p_inh = k_inh / n_inh if n_inh else 0.0
    rd = p_sys - p_inh

    fisher_primary = fisher_hypergeometric_exact(n_sys, k_sys, n_inh, k_inh)
    nw_low, nw_high = newcombe_wilson_rd_ci95(n_sys, k_sys, n_inh, k_inh)
    cz_low, cz_high = chan_zhang_unconditional_rd_ci95(n_sys, k_sys, n_inh, k_inh)

    p_inh_pure = k_inh_pure / n_inh_pure if n_inh_pure else 0.0
    fisher_pure = fisher_hypergeometric_exact(n_sys, k_sys, n_inh_pure, k_inh_pure)
    nw_pure_low, nw_pure_high = newcombe_wilson_rd_ci95(n_sys, k_sys, n_inh_pure, k_inh_pure)

    # Non-gated descriptive context: all-phase PRISMA-940
    n_sys_all = len(all_phase_systemic)
    k_sys_all = sum(1 for r in all_phase_systemic if r["has_primary_sweat_chloride"])
    n_inh_all = len(all_phase_inhaled)
    k_inh_all = sum(1 for r in all_phase_inhaled if r["has_primary_sweat_chloride"])
    fisher_all = fisher_hypergeometric_exact(n_sys_all, k_sys_all, n_inh_all, k_inh_all)

    # Pre-registered 5-part pass check (`done.md` Section 3)
    crit_margins = (n_sys >= 8) and (n_inh >= 8)
    crit_direction = (k_sys * n_inh) > (k_inh * n_sys)
    crit_fisher = fisher_primary["p_one_sided_fraction"] < Fraction(1, 100)
    crit_ci = (cz_low > 0.0) and (nw_low > 0.0)
    crit_pure_sensitivity = (
        n_inh_pure > 0
        and (k_sys * n_inh_pure) > (k_inh_pure * n_sys)
        and fisher_pure["p_one_sided_fraction"] < Fraction(1, 100)
    )
    supported = bool(crit_margins and crit_direction and crit_fisher and crit_ci and crit_pure_sensitivity)

    return {
        "benchmark_version": BENCHMARK_VERSION,
        "prediction_id": "PRED-0015",
        "evaluated_on_local": "2026-11-24",
        "done_sha256": sha256_lf_text(done_path),
        "nlm_attribution": NLM_ATTRIBUTION,
        "research_boundary": RESEARCH_BOUNDARY,
        "dataset": {
            "catalog_source_id": "clinicaltrials-gov",
            "snapshot_path": "tools/trial_atlas/snapshot.json",
            "data_timestamp": str(meta.get("data_timestamp") or meta.get("dataTimestamp") or "2026-10-09"),
            "raw_interventional_count": raw_interventional_n,
            "prisma_included_count": len(prisma_studies),
            "excluded_by_rule": excluded_counts,
            "early_phase_healthy_volunteer_excluded_count": early_hv_excluded_n,
        },
        "positive_controls": pos_controls,
        "primary_estimand_early_phase": {
            "n_sys": n_sys,
            "k_sys": k_sys,
            "p_sys": round(p_sys, 6),
            "p_sys_exact": f"{k_sys}/{n_sys}",
            "n_inh": n_inh,
            "k_inh": k_inh,
            "p_inh": round(p_inh, 6),
            "p_inh_exact": f"{k_inh}/{n_inh}",
            "risk_difference": round(rd, 6),
            "p_fisher_one_sided_exact": fisher_primary["p_one_sided_exact"],
            "p_fisher_one_sided": round(fisher_primary["p_one_sided_float"], 8),
            "p_fisher_two_sided_exact": fisher_primary["p_two_sided_exact"],
            "p_fisher_two_sided": round(fisher_primary["p_two_sided_float"], 8),
            "chan_zhang_unconditional_ci95": [round(cz_low, 6), round(cz_high, 6)],
            "newcombe_wilson_ci95": [round(nw_low, 6), round(nw_high, 6)],
            "criteria": {
                "margins_ge_8": crit_margins,
                "directional_p_sys_gt_p_inh": crit_direction,
                "fisher_p_lt_0_01": crit_fisher,
                "unconditional_and_nw_ci95_gt_0": crit_ci,
                "dual_route_sensitivity_preserved": crit_pure_sensitivity,
            },
            "outcome": "supported" if supported else "refuted",
            "ledger_standing": "open" if supported else "refuted",
            "systemic_ncts": [r["nct_id"] for r in sorted(early_systemic, key=lambda x: x["nct_id"])],
            "systemic_sweat_cl_ncts": [
                r["nct_id"] for r in sorted(early_systemic, key=lambda x: x["nct_id"]) if r["has_primary_sweat_chloride"]
            ],
            "inhaled_ncts": [r["nct_id"] for r in sorted(early_inhaled, key=lambda x: x["nct_id"])],
            "inhaled_sweat_cl_ncts": [
                r["nct_id"] for r in sorted(early_inhaled, key=lambda x: x["nct_id"]) if r["has_primary_sweat_chloride"]
            ],
            "dual_route_inhaled_ncts": [
                r["nct_id"] for r in sorted(early_inhaled, key=lambda x: x["nct_id"]) if r["dual_route_modulator_co_admin"]
            ],
        },
        "sensitivity_pure_inhaled_excluding_dual_route": {
            "n_sys": n_sys,
            "k_sys": k_sys,
            "n_inh_pure": n_inh_pure,
            "k_inh_pure": k_inh_pure,
            "p_inh_pure": round(p_inh_pure, 6),
            "p_fisher_one_sided_exact": fisher_pure["p_one_sided_exact"],
            "p_fisher_one_sided": round(fisher_pure["p_one_sided_float"], 8),
            "newcombe_wilson_ci95": [round(nw_pure_low, 6), round(nw_pure_high, 6)],
        },
        "descriptive_context_all_phases_non_gated": {
            "n_sys_all": n_sys_all,
            "k_sys_all": k_sys_all,
            "p_sys_all": round(k_sys_all / n_sys_all, 6) if n_sys_all else 0.0,
            "n_inh_all": n_inh_all,
            "k_inh_all": k_inh_all,
            "p_inh_all": round(k_inh_all / n_inh_all, 6) if n_inh_all else 0.0,
            "p_fisher_one_sided_exact": fisher_all["p_one_sided_exact"],
            "p_fisher_one_sided": round(fisher_all["p_one_sided_float"], 8),
        },
    }


def run_falsification_suite(snapshot_path: Path = DEFAULT_SNAPSHOT_JSON) -> dict[str, Any]:
    """Verify all five planted negative controls (`--falsify`)."""
    rep = evaluate_snapshot(snapshot_path)
    prim = rep["primary_estimand_early_phase"]
    n_sys, k_sys = prim["n_sys"], prim["k_sys"]
    n_inh, k_inh = prim["n_inh"], prim["k_inh"]

    # Control 1: 1,000-shuffle label-permutation null (`seed = 20261124`)
    rng = random.Random(20261124)
    outcomes_pool = [1] * (k_sys + k_inh) + [0] * ((n_sys + n_inh) - (k_sys + k_inh))
    ge_obs = 0
    null_rds: list[float] = []
    obs_rd = (k_sys / n_sys) - (k_inh / n_inh)
    for _ in range(1000):
        shuffled = list(outcomes_pool)
        rng.shuffle(shuffled)
        sim_k_sys = sum(shuffled[:n_sys])
        sim_k_inh = sum(shuffled[n_sys:])
        sim_rd = (sim_k_sys / n_sys) - (sim_k_inh / n_inh)
        null_rds.append(sim_rd)
        if sim_k_sys >= k_sys:
            ge_obs += 1
    null_rds.sort()
    perm_p = ge_obs / 1000.0
    exact_p = prim["p_fisher_one_sided"]
    null_ci_lo, null_ci_hi = null_rds[25], null_rds[974]
    if abs(perm_p - exact_p) > 0.02 or not (null_ci_lo <= 0.0 <= null_ci_hi):
        raise CompartmentalizationBenchmarkError(
            f"Falsification 1 failed: perm_p={perm_p}, exact_p={exact_p}, null_ci=[{null_ci_lo}, {null_ci_hi}]"
        )

    # Control 2: Diagnostic eligibility sweat-chloride injection blocked
    poison_diag_study = {
        "nct_id": "NCT99990001",
        "brief_title": "Phase 1 Inhaled CFTR mRNA Study (Diagnostic Sweat Chloride >= 60 mmol/L Inclusion)",
        "primary_outcomes": [
            {"measure": "Diagnostic inclusion sweat chloride >= 60 mmol/L at screening", "description": "", "time_frame": "Day 1"},
            {"measure": "Incidence of treatment-emergent adverse events", "description": "", "time_frame": "Day 28"},
        ],
    }
    diag_out = classify_primary_outcomes(poison_diag_study)
    if diag_out["has_primary_sweat_chloride"]:
        raise CompartmentalizationBenchmarkError("Falsification 2 failed: diagnostic sweat chloride leaked into primary outcome")

    # Control 3: Co-administered background modulator hijack rejected
    poison_dual = {
        "nct_id": "NCT05668741",
        "brief_title": "A Phase 1/2 Study of VX-522 Inhaled mRNA Alone and in Combination With Ivacaftor",
        "official_title": "",
        "interventions": [
            {"name": "VX-522", "type": "GENETIC"},
            {"name": "Ivacaftor", "type": "DRUG"},
        ],
        "keywords": [],
    }
    dual_cls = classify_study_regimen(poison_dual)
    if dual_cls["arm"] != "inhaled_cftr" or not dual_cls["dual_route_modulator_co_admin"]:
        raise CompartmentalizationBenchmarkError(f"Falsification 3 failed: VX-522 + ivacaftor hijacked -> {dual_cls}")

    # Control 4: Non-CFTR topical agent leak rejected
    for bad_name in ("ETD001 inhaled ENaC inhibitor", "ARO-ENaC siRNA", "BX004 nebulized bacteriophage", "Dornase alfa"):
        bad_cls = classify_study_regimen(
            {
                "nct_id": "NCT99990002",
                "brief_title": f"Phase 2 Study of {bad_name} in Cystic Fibrosis",
                "official_title": "",
                "interventions": [{"name": bad_name, "type": "DRUG"}],
                "keywords": [],
            }
        )
        if bad_cls["arm"] != "non_cftr_or_other":
            raise CompartmentalizationBenchmarkError(f"Falsification 4 failed on {bad_name}: {bad_cls}")

    # Control 5: Observational & Phase 3 / healthy-volunteer cohort leaks rejected
    obs_rec = snap.study_record(
        _to_raw_study(
            {
                "nct_id": "NCT06504589",
                "study_type": "OBSERVATIONAL",
                "conditions": ["Cystic Fibrosis"],
                "phases": [],
                "primary_outcomes": [],
            }
        )
    )
    obs_rule, _ = scope_mod.classify(obs_rec)
    if obs_rule != "X1" or is_early_phase(["PHASE2", "PHASE3"]) or is_early_phase(["PHASE3"]):
        raise CompartmentalizationBenchmarkError("Falsification 5 failed: observational or Phase 3 study leaked")
    if not is_healthy_volunteer_only("Phase 1 PK Study of VX-581 in Healthy Participants", ""):
        raise CompartmentalizationBenchmarkError("Falsification 5 failed: healthy-volunteer PK study not flagged")

    return {
        "falsify_passed": True,
        "controls_passed": 5,
        "permutation_p": round(perm_p, 4),
        "null_rd_95_interval": [round(null_ci_lo, 6), round(null_ci_hi, 6)],
        "observed_rd": round(obs_rd, 6),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT_JSON)
    ap.add_argument("--done", type=Path, default=DEFAULT_DONE_MD)
    ap.add_argument("--write-results", type=Path, default=None)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--falsify", action="store_true")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)

    if "ClinicalTrials" not in NLM_ATTRIBUTION:
        raise CompartmentalizationBenchmarkError("Missing NLM attribution")

    rep = evaluate_snapshot(args.snapshot, args.done)
    if args.write_results:
        args.write_results.write_text(json.dumps(rep, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if args.falsify:
        f_rep = run_falsification_suite(args.snapshot)
        print(
            f"PASS falsify delivery_route_compartmentalization ({f_rep['controls_passed']}/5 negative controls; "
            f"perm_p={f_rep['permutation_p']:.4f}, null_rd_95={f_rep['null_rd_95_interval']})"
        )
        return 0

    if args.check:
        f_rep = run_falsification_suite(args.snapshot)
        prim = rep["primary_estimand_early_phase"]
        print(
            f"PASS delivery_route_compartmentalization: sys={prim['k_sys']}/{prim['n_sys']} ({prim['p_sys']:.4f}) "
            f"vs inh={prim['k_inh']}/{prim['n_inh']} ({prim['p_inh']:.4f}), "
            f"p_fisher_1s={prim['p_fisher_one_sided_exact']} ({prim['p_fisher_one_sided']:.6g}), "
            f"cz_ci95={prim['chan_zhang_unconditional_ci95']}, nw_ci95={prim['newcombe_wilson_ci95']}, "
            f"outcome={prim['outcome']}, falsify={f_rep['controls_passed']}/5"
        )
        return 0

    if args.json:
        print(json.dumps(rep, indent=2, ensure_ascii=False))
    else:
        prim = rep["primary_estimand_early_phase"]
        print(
            f"PRED-0015 ({prim['outcome']}): systemic={prim['k_sys']}/{prim['n_sys']} "
            f"vs inhaled={prim['k_inh']}/{prim['n_inh']}, p_fisher_1s={prim['p_fisher_one_sided_exact']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
