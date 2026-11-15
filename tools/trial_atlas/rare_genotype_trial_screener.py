#!/usr/bin/env python3
"""Deterministic Trial Endpoint Match & Eligibility Screener for Rare/Class-I CF Genotypes (T-0920).

Screens PRISMA-included interventional Cystic Fibrosis clinical trials in
`tools/trial_atlas/snapshot.json` for rare/Class I stop-codon (`G542X`, `W1282X`,
`R553X`, `R1162X`), frameshift, splicing (`3849+10kbC->T`, `2789+5G->A`), and
genotype-agnostic modalities (mRNA/LNP, gene therapy/editing, readthrough/tRNA,
phage/anti-infective) alongside their stated clinical endpoints (`sweat_chloride`,
`ppfev1`, `lci`, `safety_tolerability`, `pk_pd`, `microbiology_cfu`).

Safety & Grounding Rules (`cf-ai-safe-use`, `cf-visit-prep`, `cf-research-context`):
- 100% pure Python 3.10+ standard library, zero network calls, zero LLM calls.
- Refuses personal health information (PHI: dates of birth, MRNs, emails, long IDs)
  via `care_team_trials.validate_query_inputs`.
- Never makes personal clinical eligibility promises ("you are eligible",
  "qualifies for", "is a candidate"). Every row carries `[T1: ClinicalTrials.gov]`
  provenance and quotes verbatim registry fields (`nct_id`, `brief_title`,
  `overall_status`, `phases`, `modalities`, and stated primary outcome measures),
  reporting `"not stated"` where the registry record is silent.
- Bilingual (`en` and `es-MX`) clinic-visit discussion questions for care teams.
- Deterministic `--falsify` mode runs poisoned negative controls (PHI input,
  personal eligibility claim injection, and ungrounded NCT ID injection).
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from care_team_trials import (  # noqa: E402
    DEFAULT_SNAPSHOT_JSON,
    NOT_STATED,
    RefusedInputError,
    load_prisma_included_studies,
    validate_query_inputs,
)

DOSE_SHORTHAND_RE = re.compile(r"\b\d+(?:\.\d+)?\s*(?:mg|mcg|ug|g|ml|iu)\b", re.I)


def sanitize_registry_text(raw: Any) -> str:
    """Clean registry text and strip medication dose shorthand per cf-ai-safe-use."""
    text = str(raw or "").strip()
    if not text:
        return NOT_STATED
    text = DOSE_SHORTHAND_RE.sub("[dose omitted]", text)
    return re.sub(r"\s+", " ", text).strip()

VARIANT_CLASS_RULES: dict[str, tuple[str, ...]] = {
    "class_i_stop_nonsense": (
        "nonsense",
        "stop codon",
        "premature termination",
        "ptc",
        "readthrough",
        "g542x",
        "w1282x",
        "r553x",
        "r1162x",
        "y122x",
        "elx-02",
        "ataluren",
        "trna",
        "class i",
    ),
    "splice_modulation": (
        "splic",
        "3849+10kb",
        "2789+5g",
        "antisense",
        "oligonucleotide",
        "aso",
        "spl84",
    ),
    "genotype_agnostic_mrna_gene": (
        "mrna",
        "gene therapy",
        "gene editing",
        "genetic",
        "lentivir",
        "aav",
        "lipid nanoparticle",
        "lnp",
        "vx-522",
        "arct-032",
        "4d-710",
        "bi 3720931",
        "rct2100",
        "inhalation of cftr",
        "modulator-ineligible",
        "not eligible for cftr modulator",
    ),
    "phage_anti_infective": (
        "phage",
        "bacteriophage",
        "pseudomonas",
        "anti-infective",
        "inhaled antibiotic",
        "biofilm",
    ),
}

ENDPOINT_CATEGORY_RULES: dict[str, tuple[str, ...]] = {
    "sweat_chloride": ("sweat chloride", "sweat cl"),
    "ppfev1": ("fev1", "ppfev1", "forced expiratory volume", "spirometry"),
    "lci": ("lung clearance index", "lci", "lci2.5", "multiple breath washout"),
    "safety_tolerability": ("adverse event", "tea", "teae", "safety", "tolerability", "dose-limiting"),
    "pk_pd": ("pharmacokinetic", "pk", "cmax", "auc", "mrna expression", "cftr protein"),
    "microbiology_cfu": ("colony forming", "cfu", "sputum density", "pseudomonas", "bacterial load"),
}

FORBIDDEN_ELIGIBILITY_PATTERNS: tuple[tuple[str, re.Pattern[str]], ...] = (
    ("personal_eligibility", re.compile(r"\b(?:you\s+are\s+eligible|you\s+may\s+be\s+eligible|you\s+qualify|qualifies\s+for\s+this\s+trial|is\s+a\s+candidate\s+for)\b", re.I)),
    ("clinical_recommendation", re.compile(r"\b(?:you\s+should\s+enroll|we\s+recommend\s+enrolling|guaranteed\s+benefit|curative\s+guarantee)\b", re.I)),
)


class ScreenerVerificationError(ValueError):
    """Raised when an output record violates safety or registry provenance rules."""


def _study_corpus(study: dict[str, Any]) -> str:
    parts: list[str] = [
        str(study.get("brief_title") or ""),
        " ".join(str(m) for m in (study.get("modalities") or [])),
        " ".join(str(i) for i in (study.get("interventions") or [])),
        " ".join(
            f"{o.get('measure', '')} {o.get('time_frame', '')}"
            for o in (study.get("primary_outcomes") or [])
            if isinstance(o, dict)
        ),
    ]
    return " ".join(parts).lower()


def classify_variant_classes(study: dict[str, Any]) -> list[str]:
    """Infer stated trial relevance categories from registry title, modalities, and interventions."""
    corpus = _study_corpus(study)
    matched: list[str] = []
    for cls_name, keywords in VARIANT_CLASS_RULES.items():
        if any(kw in corpus for kw in keywords):
            matched.append(cls_name)
    for mod in study.get("modalities") or []:
        mod_l = str(mod).lower()
        if any(k in mod_l for k in ("mrna", "gene", "genetic", "rna", "editing")) and "genotype_agnostic_mrna_gene" not in matched:
            matched.append("genotype_agnostic_mrna_gene")
        if "phage" in mod_l and "phage_anti_infective" not in matched:
            matched.append("phage_anti_infective")
    return sorted(matched)


def classify_endpoint_categories(study: dict[str, Any]) -> list[str]:
    """Classify stated primary outcome measures into standardized clinical endpoint categories."""
    outcomes = study.get("primary_outcomes") or []
    outcome_text = " ".join(
        str(o.get("measure") or "") for o in outcomes if isinstance(o, dict)
    ).lower()
    matched: list[str] = []
    for ep_name, keywords in ENDPOINT_CATEGORY_RULES.items():
        if any(kw in outcome_text for kw in keywords):
            matched.append(ep_name)
    return sorted(matched)


def screen_rare_genotype_trials(
    studies: list[dict[str, Any]],
    variant_class: str | None = None,
    variant_symbol: str | None = None,
    status: str | None = "RECRUITING",
    endpoint: str | None = None,
    limit: int = 15,
) -> list[dict[str, Any]]:
    """Deterministically screen PRISMA-included trials by rare/Class-I criteria and endpoints."""
    validate_query_inputs(variant_class, variant_symbol, status, endpoint)

    vc_norm = variant_class.strip().lower() if variant_class and variant_class.strip().lower() != "all" else None
    if vc_norm and vc_norm not in VARIANT_CLASS_RULES:
        raise ValueError(
            f"Unknown --variant-class '{variant_class}'. Allowed: {sorted(VARIANT_CLASS_RULES)} or 'all'."
        )

    ep_norm = endpoint.strip().lower() if endpoint and endpoint.strip().lower() != "all" else None
    if ep_norm and ep_norm not in ENDPOINT_CATEGORY_RULES:
        raise ValueError(
            f"Unknown --endpoint '{endpoint}'. Allowed: {sorted(ENDPOINT_CATEGORY_RULES)} or 'all'."
        )

    var_norm = variant_symbol.strip().lower() if variant_symbol and variant_symbol.strip() else None
    status_norm = status.strip().upper() if status and status.strip().upper() != "ALL" else None

    rows: list[dict[str, Any]] = []
    for s in studies:
        s_status = str(s.get("overall_status") or NOT_STATED).upper()
        if status_norm and s_status != status_norm:
            continue

        v_classes = classify_variant_classes(s)
        if vc_norm and vc_norm not in v_classes:
            continue
        if not vc_norm and not v_classes:
            continue

        ep_cats = classify_endpoint_categories(s)
        if ep_norm and ep_norm not in ep_cats:
            continue

        corpus = _study_corpus(s)
        if var_norm and var_norm not in corpus:
            # Allow matching genotype-agnostic genetic/mRNA trials when screening a Class I/rare variant
            if "genotype_agnostic_mrna_gene" not in v_classes and "class_i_stop_nonsense" not in v_classes:
                continue

        outcomes = s.get("primary_outcomes") or []
        first_outcome = outcomes[0] if outcomes and isinstance(outcomes[0], dict) else {}
        primary_measure = sanitize_registry_text(first_outcome.get("measure")) if first_outcome else NOT_STATED
        primary_timeframe = sanitize_registry_text(first_outcome.get("time_frame")) if first_outcome else NOT_STATED

        rows.append(
            {
                "nct_id": str(s.get("nct_id") or NOT_STATED),
                "provenance_tier": "T1",
                "provenance_source": "ClinicalTrials.gov API v2 snapshot (PRISMA-included)",
                "brief_title": sanitize_registry_text(s.get("brief_title")),
                "overall_status": str(s.get("overall_status") or NOT_STATED),
                "phases": list(s.get("phases") or [NOT_STATED]),
                "modalities": list(s.get("modalities") or [NOT_STATED]),
                "matched_variant_classes": v_classes,
                "matched_endpoint_categories": ep_cats if ep_cats else [NOT_STATED],
                "stated_primary_outcome_measure": primary_measure,
                "stated_primary_outcome_time_frame": primary_timeframe,
                "explicit_variant_in_title_or_intervention": (
                    variant_symbol.strip()
                    if (var_norm and var_norm in corpus and variant_symbol)
                    else NOT_STATED
                ),
            }
        )
        if len(rows) >= max(1, limit):
            break

    return rows


def build_care_team_packet(
    rows: list[dict[str, Any]],
    known_ncts: set[str],
    lang: str = "en",
    filters: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a verified bilingual care-team discussion packet with mandatory safety disclaimers."""
    is_es = lang.strip().lower().startswith("es")
    if is_es:
        boundary = (
            "Solo para preparacion informativa de consulta (Tier T1: ClinicalTrials.gov). "
            "Esta herramienta NO evalua elegibilidad personal ni brinda consejo medico."
        )
        questions = [
            "¿Evalua este estudio variantes de Clase I / codones de paro o terapias independientes del genotipo (ARNm/genetica)?",
            "¿Cuales son los criterios de inclusion/exclusion exactos y que mediciones primarias (FEV1, cloruro en sudor, LCI) se requieren?",
            "¿Existe algun centro del estudio cerca o apoyo de traslado a traves de la fundacion de FQ?",
        ]
        urgent_line = (
            "Si presenta dificultad respiratoria aguda, hemoptisis o fiebre alta, comuniquese "
            "de inmediato con su equipo de atencion de fibrosis quistica o servicios de urgencias."
        )
    else:
        boundary = (
            "Informational clinic-visit preparation only (Tier T1: ClinicalTrials.gov). "
            "This tool does NOT determine personal trial eligibility or provide medical advice."
        )
        questions = [
            "Does this registered study enroll people with Class I / nonsense / splice variants or evaluate a genotype-agnostic therapy?",
            "What are the exact protocol inclusion/exclusion criteria and stated primary endpoints (ppFEV1, sweat chloride, LCI, or safety)?",
            "Are any study sites open near us, and how does trial participation interact with current airway clearance or standard care?",
        ]
        urgent_line = (
            "If you have urgent symptoms such as sudden shortness of breath, coughing up blood, "
            "or high fever, contact your CF care team or emergency services right away."
        )

    packet = {
        "schema_version": 1,
        "tool": "rare_genotype_trial_screener",
        "lang": "es-MX" if is_es else "en",
        "epistemic_tier": "T1",
        "clinical_boundary": boundary,
        "urgent_symptom_notice": urgent_line,
        "filters": filters or {},
        "match_count": len(rows),
        "trials": rows,
        "care_team_questions": questions,
    }
    verify_packet_safety_and_provenance(packet, known_ncts)
    return packet


def verify_packet_safety_and_provenance(packet: dict[str, Any], known_ncts: set[str]) -> None:
    """Deterministic verifier enforcing zero personal eligibility promises and strict NCT provenance."""
    serialized = json.dumps(packet, ensure_ascii=False)
    for rule_name, pat in FORBIDDEN_ELIGIBILITY_PATTERNS:
        if pat.search(serialized):
            raise ScreenerVerificationError(
                f"Safety violation ({rule_name}): output contains forbidden personal eligibility or clinical promise phrasing."
            )
    if not packet.get("clinical_boundary") or "not" not in str(packet.get("clinical_boundary", "")).lower() and "no " not in str(packet.get("clinical_boundary", "")).lower():
        raise ScreenerVerificationError("Safety violation: missing explicit non-medical-advice clinical_boundary.")

    for row in packet.get("trials", []):
        nct = str(row.get("nct_id") or "")
        if not re.fullmatch(r"NCT\d{8}", nct):
            raise ScreenerVerificationError(f"Provenance violation: invalid NCT ID '{nct}'.")
        if known_ncts and nct not in known_ncts:
            raise ScreenerVerificationError(
                f"Provenance violation: NCT ID '{nct}' is absent from the verified PRISMA snapshot."
            )
        if row.get("provenance_tier") != "T1":
            raise ScreenerVerificationError(
                f"Provenance violation: trial '{nct}' missing provenance_tier='T1'."
            )


def run_falsification_suite(snapshot_path: Path = DEFAULT_SNAPSHOT_JSON) -> dict[str, Any]:
    """Execute deterministic negative controls (--falsify) and verify every poisoned case fails closed."""
    _, studies = load_prisma_included_studies(snapshot_path)
    known_ncts = {str(s.get("nct_id")) for s in studies if s.get("nct_id")}
    controls: list[dict[str, Any]] = []

    # Negative Control 1: PHI query input must raise RefusedInputError
    try:
        validate_query_inputs("Patient DOB 04/12/2014 MRN 99481234 G542X")
        controls.append({"control": "phi_input_refusal", "passed": False, "detail": "PHI input was not refused"})
    except RefusedInputError as exc:
        controls.append({"control": "phi_input_refusal", "passed": True, "detail": str(exc)})

    # Negative Control 2: Poisoned personal eligibility claim in packet must raise ScreenerVerificationError
    poisoned_eligibility = {
        "clinical_boundary": "Informational only, not medical advice.",
        "trials": [],
        "care_team_questions": ["Based on G542X, you are eligible and we recommend enrolling immediately."],
    }
    try:
        verify_packet_safety_and_provenance(poisoned_eligibility, known_ncts)
        controls.append({"control": "personal_eligibility_promise_guard", "passed": False, "detail": "Eligibility promise was not caught"})
    except ScreenerVerificationError as exc:
        controls.append({"control": "personal_eligibility_promise_guard", "passed": True, "detail": str(exc)})

    # Negative Control 3: Hallucinated NCT ID absent from snapshot must raise ScreenerVerificationError
    poisoned_nct = {
        "clinical_boundary": "Informational only, not medical advice.",
        "trials": [
            {
                "nct_id": "NCT99999999",
                "provenance_tier": "T1",
                "brief_title": "Hallucinated Unregistered Trial",
            }
        ],
    }
    try:
        verify_packet_safety_and_provenance(poisoned_nct, known_ncts)
        controls.append({"control": "hallucinated_nct_provenance_guard", "passed": False, "detail": "Hallucinated NCT99999999 was not rejected"})
    except ScreenerVerificationError as exc:
        controls.append({"control": "hallucinated_nct_provenance_guard", "passed": True, "detail": str(exc)})

    all_passed = all(c["passed"] for c in controls)
    return {
        "tool": "rare_genotype_trial_screener",
        "falsify_passed": all_passed,
        "negative_controls": controls,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT_JSON, help="Path to trial_atlas/snapshot.json")
    parser.add_argument("--variant-class", type=str, default="all", help="Variant/modality class filter (e.g. genotype_agnostic_mrna_gene, class_i_stop_nonsense, splice_modulation, phage_anti_infective, all)")
    parser.add_argument("--variant", type=str, default="", help="Optional variant symbol (e.g. G542X, W1282X, 3849+10kbC->T)")
    parser.add_argument("--status", type=str, default="RECRUITING", help="Registry recruitment status (default: RECRUITING, or ALL)")
    parser.add_argument("--endpoint", type=str, default="all", help="Endpoint category filter (sweat_chloride, ppfev1, lci, safety_tolerability, pk_pd, microbiology_cfu, all)")
    parser.add_argument("--lang", type=str, default="en", choices=("en", "es-MX"), help="Output language (en or es-MX)")
    parser.add_argument("--limit", type=int, default=15, help="Max trials to return")
    parser.add_argument("--json", action="store_true", help="Emit JSON packet")
    parser.add_argument("--check", action="store_true", help="Run fast deterministic self-check")
    parser.add_argument("--falsify", action="store_true", help="Run negative-control falsification suite")
    args = parser.parse_args(argv)

    if args.falsify:
        rep = run_falsification_suite(args.snapshot)
        if args.json:
            print(json.dumps(rep, indent=2))
        else:
            print(f"PASS falsify rare_genotype_trial_screener ({len(rep['negative_controls'])}/3 negative controls)" if rep["falsify_passed"] else "FAIL falsify rare_genotype_trial_screener")
        return 0 if rep["falsify_passed"] else 1

    _, studies = load_prisma_included_studies(args.snapshot)
    known_ncts = {str(s.get("nct_id")) for s in studies if s.get("nct_id")}
    rows = screen_rare_genotype_trials(
        studies,
        variant_class=args.variant_class,
        variant_symbol=args.variant or None,
        status=args.status,
        endpoint=args.endpoint,
        limit=args.limit,
    )
    packet = build_care_team_packet(
        rows,
        known_ncts=known_ncts,
        lang=args.lang,
        filters={
            "variant_class": args.variant_class,
            "variant": args.variant or NOT_STATED,
            "status": args.status,
            "endpoint": args.endpoint,
        },
    )

    if args.check:
        f_rep = run_falsification_suite(args.snapshot)
        if not f_rep["falsify_passed"] or packet["match_count"] < 1:
            print("FAIL rare_genotype_trial_screener self-check", file=sys.stderr)
            return 1
        print(
            f"PASS rare_genotype_trial_screener: prisma_studies={len(studies)} "
            f"matched={packet['match_count']} falsify_controls={len(f_rep['negative_controls'])}/3"
        )
        return 0

    if args.json:
        print(json.dumps(packet, indent=2, ensure_ascii=False))
    else:
        print(f"=== {packet['tool']} ({packet['lang']}, {packet['epistemic_tier']}) ===")
        print(f"Boundary: {packet['clinical_boundary']}")
        print(f"Matched Trials: {packet['match_count']}")
        for r in packet["trials"]:
            print(
                f"  - [{r['nct_id']}] ({r['overall_status']}, {','.join(r['phases'])}) "
                f"{r['brief_title']} | classes={','.join(r['matched_variant_classes'])} "
                f"| primary_endpoint={r['stated_primary_outcome_measure']}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
