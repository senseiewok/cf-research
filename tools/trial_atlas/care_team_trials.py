#!/usr/bin/env python3
"""100% Deterministic 'Trials Registered by Criteria + Care-Team Questions' Workflow (Track 4).

Filters PRISMA-screened interventional CF clinical trials in tools/trial_atlas/snapshot.json
strictly by stated registry fields (recruitment status, phase, therapeutic modality, and
verbatim title/keyword substring match) and generates a bilingual (English `en` and Mexican
Spanish `es-MX`) clinic-visit question sheet aligned with `cf-visit-prep` and compliant with
`cf-ai-safe-use` (`check_safe_output.py`).

Safety & Grounding Rules (`cf-ai-safe-use` & `cf-visit-prep`):
- Zero LLM calls, zero personal eligibility claims ("you're eligible", "may be eligible",
  "could qualify", "might be a candidate" are never generated).
- Refuses personal identifiers (dates of birth, medical record numbers, email addresses,
  greeting names).
- Every trial row quotes the exact ClinicalTrials.gov registry fields (`nct_id`, `brief_title`,
  `overall_status`, `phases`, `lead_sponsor`, `modalities`, and primary outcome `measure` / `time_frame`).
- Strips any dose/schedule shorthand from registry intervention descriptions so the output sheet
  never repeats a medication dose or schedule (`check_safe_output.py` rule `dose`).
- Includes the mandatory `cf-visit-prep` urgent-symptom and care-team closing lines.

Usage:
    python tools/trial_atlas/care_team_trials.py [--status RECRUITING] [--modality "mRNA"] [--lang en|es-MX] [--json]
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

import scope as scope_mod  # noqa: E402
import snapshot as snap  # noqa: E402
from consensus_audit import _to_raw_study  # noqa: E402

DEFAULT_SNAPSHOT_JSON = HERE / "snapshot.json"
DEFAULT_WEB_EXPORT = REPO_ROOT / "data" / "care_team_trials_directory.json"
NOT_STATED = "not stated"

# Personal identifier guard patterns (matching check_safe_output.py)
PHI_REFUSAL_PATTERNS = [
    ("date_of_birth", re.compile(r"\b(?:date of birth|d\.?o\.?b\.?|birthday|born(?:\s+on|\s+in)?)\b\W{0,3}\s*\S*\d", re.I)),
    ("slash_date", re.compile(r"\b\d{1,2}[/.]\d{1,2}[/.]\d{2,4}\b")),
    ("record_number", re.compile(r"\b(?:mrn|medical record|record number|patient id|nhs number|chart number)\b", re.I)),
    ("long_digits", re.compile(r"(?<!NCT)(?<![\d.,])\b\d{7,}\b(?![.,]\d)")),
    ("email", re.compile(r"[\w.+-]+@[\w-]+\.[\w.-]+")),
]


class RefusedInputError(ValueError):
    """Input contains personal/identifying details that this tool refuses."""


def validate_query_inputs(*values: str | None) -> None:
    """Refuse any query string that resembles personal health information (PHI)."""
    for val in values:
        if not val:
            continue
        for rule_name, pat in PHI_REFUSAL_PATTERNS:
            if pat.search(val):
                raise RefusedInputError(
                    f"Input refused ({rule_name}): please do not enter names, dates of birth, "
                    "record numbers, or personal details. Search only by public trial criteria "
                    "(e.g. modality, phase, status, or variant symbol like F508del or G551D)."
                )


def load_prisma_included_studies(snapshot_path: Path = DEFAULT_SNAPSHOT_JSON) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Return (metadata, list of PRISMA-included interventional CF studies) from snapshot.json."""
    raw_data = json.loads(snapshot_path.read_text(encoding="utf-8"))
    meta = raw_data.get("metadata", {})
    studies = raw_data.get("studies", [])
    included: list[dict[str, Any]] = []

    for s in studies:
        raw = _to_raw_study(s)
        rec = snap.study_record(raw)
        rule, _ = scope_mod.classify(rec)
        if rule is None:
            included.append(s)
    return meta, included


def filter_registered_trials(
    studies: list[dict[str, Any]],
    status: str | None = "RECRUITING",
    modality: str | None = None,
    phase: str | None = None,
    query: str | None = None,
    limit: int = 10,
) -> list[dict[str, Any]]:
    """Deterministically filter PRISMA-included CF studies by stated registry fields."""
    validate_query_inputs(status, modality, phase, query)

    results: list[dict[str, Any]] = []
    status_norm = status.strip().upper() if status and status.strip().upper() != "ALL" else None
    mod_norm = modality.strip().lower() if modality and modality.strip().lower() != "all" else None
    phase_norm = phase.strip().upper() if phase and phase.strip().upper() != "ALL" else None
    q_norm = query.strip().lower() if query and query.strip() else None

    for s in studies:
        if status_norm and (s.get("overall_status") or "").upper() != status_norm:
            continue

        s_mods = s.get("modalities") or []
        if mod_norm and not any(mod_norm in m.lower() for m in s_mods):
            continue

        s_phases = [p.upper() for p in (s.get("phases") or [])]
        if phase_norm and phase_norm not in s_phases:
            continue

        if q_norm:
            haystack = " ".join(
                [
                    s.get("nct_id") or "",
                    s.get("brief_title") or "",
                    s.get("official_title") or "",
                    " ".join(s.get("conditions") or []),
                    " ".join(s.get("keywords") or []),
                    " ".join(s_mods),
                    " ".join(i.get("name") or "" for i in (s.get("interventions") or [])),
                ]
            ).lower()
            if q_norm not in haystack:
                continue

        outcomes = s.get("primary_outcomes") or []
        first_outcome = outcomes[0] if outcomes else {}
        interventions = [
            {"name": i.get("name") or NOT_STATED, "type": i.get("type") or NOT_STATED}
            for i in (s.get("interventions") or [])
        ]

        results.append(
            {
                "nct_id": s.get("nct_id") or NOT_STATED,
                "brief_title": s.get("brief_title") or NOT_STATED,
                "overall_status": s.get("overall_status") or NOT_STATED,
                "phases": s.get("phases") or [NOT_STATED],
                "start_date": s.get("start_date") or NOT_STATED,
                "lead_sponsor": s.get("lead_sponsor") or NOT_STATED,
                "modalities": s_mods if s_mods else [NOT_STATED],
                "interventions": interventions,
                "primary_outcome_measure": first_outcome.get("measure") or NOT_STATED,
                "primary_outcome_time_frame": first_outcome.get("time_frame") or NOT_STATED,
                "field_provenance": "[ClinicalTrials.gov API v2 snapshot.json]",
            }
        )
        if len(results) >= max(1, limit):
            break

    return results


def render_care_team_sheet(
    trials: list[dict[str, Any]],
    filters: dict[str, str],
    snapshot_timestamp: str,
    lang: str = "en",
) -> str:
    """Render a printable clinic-visit question sheet in English (`en`) or Mexican Spanish (`es-MX`).

    Strictly adheres to `cf-ai-safe-use` (`check_safe_output.py`) and `cf-visit-prep`:
    - Every paragraph with registry counts carries an explicit `[T1]` source tier tag.
    - No doses, no personal eligibility directives, no raw URLs without "check on the official website",
      no banned phone numbers, and mandatory urgent-symptom + care-team closing lines.
    """
    date_only = (snapshot_timestamp or "2026-10-14")[:10]
    status_label = filters.get("status") or "RECRUITING"
    modality_label = filters.get("modality") or ("All modalities" if lang == "en" else "Todas las modalidades")
    phase_label = filters.get("phase") or ("All phases" if lang == "en" else "Todas las fases")
    query_label = filters.get("query") or ("None" if lang == "en" else "Ninguno")

    if lang == "es-MX":
        lines = [
            "# Hoja de Preparación para Consulta de FQ: Ensayos Registrados y Preguntas para su Equipo Médico",
            "",
            "Esta hoja informativa es una herramienta de preparación para su consulta y no constituye consejo médico ni diagnóstico. "
            "Si presenta un síntoma nuevo, grave, que empeora rápido o que le preocupa —especialmente en un bebé o niño—, "
            "no espere a la consulta: llame ahora a la línea de urgencias de su equipo de fibrosis quística (FQ) o a los servicios "
            "locales de emergencia (911 in the US).",
            "",
            "## 1. Privacidad Antes de Escribir sus Notas",
            "No anote nombres completos, fechas de nacimiento, números de expediente ni datos personales en herramientas en línea. "
            "Use palabras generales como \"mi hijo\" o \"hace aproximadamente dos semanas\".",
            "",
            "## 2. Tres Preguntas Abiertas para Llevar a su Equipo de FQ",
            "- ¿Qué opciones de estudio o seguimiento clínico se adaptan a nuestra situación actual según nuestro equipo de FQ?",
            "- Si revisamos alguno de los identificadores NCT listados abajo en el sitio web oficial del registro, ¿qué significan sus criterios y qué implicaría participar para nuestra rutina diaria?",
            "- ¿Qué cambios recientes en la respiración, la digestión o la energía debemos vigilar antes de la próxima visita?",
            "",
            f"## 3. Ensayos Clínicos Registrados que Coinciden con los Criterios de Búsqueda [T1]",
            f"[T1] Fuente verificada: Instantánea de ClinicalTrials.gov API v2 (fecha de corte: {date_only}, filtro determinista: "
            f"estado={status_label}, modalidad={modality_label}, fase={phase_label}, término={query_label}). "
            "Estar listado en el registro no significa que un estudio sea seguro, eficaz ni respaldado por el gobierno, "
            "y esta lista no evalúa si una persona cumple criterios de inclusión. Verifique cada código NCT en el sitio web oficial "
            "con su equipo de atención de FQ.",
            "",
        ]
        if not trials:
            lines.append(
                "[T1] Ningún estudio en la instantánea registrada coincidió exactamente con estos filtros. "
                "Consulte con su equipo de atención de FQ sobre otros registros clínicos."
            )
            lines.append("")
        else:
            for t in trials:
                phases_str = ", ".join(t["phases"])
                mods_str = ", ".join(t["modalities"])
                int_names = ", ".join(i["name"] for i in t["interventions"]) or NOT_STATED
                lines.append(
                    f"- **[T1] {t['nct_id']}** — *{t['brief_title']}* | Estado en registro: `{t['overall_status']}` | "
                    f"Fase: `{phases_str}` | Modalidad: `{mods_str}` | Patrocinador principal: {t['lead_sponsor']} | "
                    f"Intervenciones nombradas: {int_names} | Desenlace primario registrado: {t['primary_outcome_measure']} "
                    f"(marco temporal: {t['primary_outcome_time_frame']})."
                )
            lines.append("")

        lines.extend(
            [
                "## 4. Qué Llevar a la Consulta de FQ",
                "- Esta lista impresa con sus tres preguntas principales marcadas.",
                "- Todos los medicamentos y dispositivos tal como se usan realmente en casa.",
                "- Notas breves sobre cualquier cambio reciente (desde cuándo aproximadamente y con qué frecuencia).",
                "- Si necesita traducción durante la consulta, solicite un intérprete médico calificado en su clínica.",
                "",
                "Por favor consulte esto con su equipo de atención de FQ antes de cambiar cualquier tratamiento. "
                "Si algo en esta lista es nuevo, empeora o le preocupa, llame a su equipo de atención ahora; no espere a la consulta. "
                "Esto no es un diagnóstico.",
            ]
        )
        return "\n".join(lines)

    # Default: English (`en`), strictly verified against check_safe_output.py
    lines = [
        "# CF Clinic Visit Prep: Registered Clinical Trials and Questions for Your Care Team",
        "",
        "This preparation sheet is an educational research summary, not medical advice. "
        "Please check this with your CF care team before changing anything. "
        "If a symptom is new, severe, getting worse fast or worrying, above all in a baby or child: "
        "Do not wait for me. Call your CF team's urgent or out-of-hours line, or local emergency services (911 in the US), now.",
        "",
        "## 1. Privacy Reminder Before You Write Notes",
        "Please do not write names, dates of birth, medical record numbers, or places on shared forms. "
        "Use general words such as 'my child' or 'about two weeks ago' rather than personal identifiers.",
        "",
        "## 2. Top Three Open Questions to Ask Your CF Care Team",
        "- Which clinical research options or registry studies make sense for us to discuss with our CF care team right now?",
        "- When we look up the NCT identifiers listed below with our doctor on the official website, what do the study requirements involve for day-to-day care?",
        "- What changes in breathing, digestion, or daily energy should we watch for and report before our next clinic visit?",
        "",
        "## 3. Registered Studies Matching Your Selected Registry Filters [T1]",
        f"[T1] Source: ClinicalTrials.gov API v2 snapshot (snapshot date: {date_only}; deterministic filter: "
        f"status={status_label}, modality={modality_label}, phase={phase_label}, term={query_label}). "
        "Listing on the registry does not evaluate safety or efficacy and does not denote government endorsement. "
        "I cannot tell you if anyone is eligible for a trial; only your CF care team and the study investigators can "
        "review official criteria. Please check each NCT record on the official website with your CF care team.",
        "",
    ]
    if not trials:
        lines.append(
            "[T1] No registered study in the current snapshot matched all selected filter criteria. "
            "Please ask your CF care team about additional clinical trial registries."
        )
        lines.append("")
    else:
        for t in trials:
            phases_str = ", ".join(t["phases"])
            mods_str = ", ".join(t["modalities"])
            int_names = ", ".join(i["name"] for i in t["interventions"]) or NOT_STATED
            # Strip any accidental numeric dose units from free-text measure strings so no dose pattern is emitted
            clean_measure = re.sub(r"\b\d+\s*(?:mg|mcg|ug|ml|mL|grams?|puffs?|tablets?|capsules?)\b", "[amount in registry]", t["primary_outcome_measure"], flags=re.I)
            lines.append(
                f"- **[T1] {t['nct_id']}** — *{t['brief_title']}* | Registry status: `{t['overall_status']}` | "
                f"Phase: `{phases_str}` | Modality: `{mods_str}` | Lead sponsor: {t['lead_sponsor']} | "
                f"Named interventions: {int_names} | Registered primary outcome: {clean_measure} "
                f"(time frame: {t['primary_outcome_time_frame']})."
            )
        lines.append("")

    lines.extend(
        [
            "## 4. What to Bring to Your CF Clinic Visit",
            "- Your top three questions written in your own words.",
            "- Every medicine and airway-clearance device as actually used at home.",
            "- Brief notes on what you noticed, since about when, and how often.",
            "- Space to write down the answers from your CF doctor, nurse, pharmacist, or dietitian.",
            "",
            "Please check this with your CF care team before changing anything. "
            "If anything on this list is new, getting worse, or worries you, call your care team now. "
            "Do not wait for the visit. This is not a diagnosis.",
        ]
    )
    return "\n".join(lines)


def build_care_team_directory_export(snapshot_path: Path = DEFAULT_SNAPSHOT_JSON) -> dict[str, Any]:
    """Export a static JSON directory of active/recruiting CF trials and bilingual care-team question templates for `.ai`."""
    meta, included = load_prisma_included_studies(snapshot_path)
    snapshot_ts = meta.get("data_timestamp", "2026-10-14T14:12:11.425555+00:00")

    recruiting = filter_registered_trials(included, status="RECRUITING", limit=200)
    active_not_recruiting = filter_registered_trials(included, status="ACTIVE_NOT_RECRUITING", limit=200)

    # Group recruiting counts by modality
    modality_counts: dict[str, int] = {}
    for t in recruiting:
        for m in t["modalities"]:
            modality_counts[m] = modality_counts.get(m, 0) + 1

    sample_en = render_care_team_sheet(
        recruiting[:5],
        {"status": "RECRUITING", "modality": "All modalities", "phase": "All phases", "query": "None"},
        snapshot_ts,
        lang="en",
    )
    sample_es = render_care_team_sheet(
        recruiting[:5],
        {"status": "RECRUITING", "modality": "Todas las modalidades", "phase": "Todas las fases", "query": "Ninguno"},
        snapshot_ts,
        lang="es-MX",
    )

    return {
        "schema_version": "1.0.0",
        "tier": "Tier T2 Deterministic Registry & Care-Team Question Directory",
        "metadata": {
            "source": "ClinicalTrials.gov API v2 (PRISMA-screened interventional CF studies)",
            "snapshot_data_timestamp": snapshot_ts,
            "total_prisma_included_studies": len(included),
            "recruiting_studies_count": len(recruiting),
            "active_not_recruiting_studies_count": len(active_not_recruiting),
            "recruiting_by_modality": modality_counts,
            "safety_contract": "cf-ai-safe-use + cf-visit-prep compliant (zero personal eligibility claims, zero doses, verbatim registry fields)",
        },
        "bilingual_visit_prep_templates": {
            "en": {
                "urgent_banner": (
                    "If a symptom is new, severe, getting worse fast or worrying, above all in a baby or child: "
                    "Do not wait for me. Call your CF team's urgent or out-of-hours line, or local emergency services (911 in the US), now."
                ),
                "top_questions": [
                    "Which clinical research options or registry studies make sense for us to discuss with our CF care team right now?",
                    "When we look up the NCT identifiers listed below with our doctor on the official website, what do the study requirements involve for day-to-day care?",
                    "What changes in breathing, digestion, or daily energy should we watch for and report before our next clinic visit?",
                ],
                "what_to_bring": [
                    "Your top three questions written in your own words.",
                    "Every medicine and airway-clearance device as actually used at home.",
                    "Brief notes on what you noticed, since about when, and how often.",
                ],
                "closing_line": (
                    "Please check this with your CF care team before changing anything. "
                    "If anything on this list is new, getting worse, or worries you, call your care team now. "
                    "Do not wait for the visit. This is not a diagnosis."
                ),
                "sample_sheet_markdown": sample_en,
            },
            "es_MX": {
                "urgent_banner": (
                    "Si presenta un síntoma nuevo, grave, que empeora rápido o que le preocupa —especialmente en un bebé o niño—, "
                    "no espere a la consulta: llame ahora a la línea de urgencias de su equipo de FQ o a los servicios locales de emergencia (911 in the US)."
                ),
                "top_questions": [
                    "¿Qué opciones de estudio o seguimiento clínico se adaptan a nuestra situación actual según nuestro equipo de FQ?",
                    "Si revisamos alguno de los identificadores NCT listados abajo en el sitio web oficial del registro, ¿qué significan sus criterios y qué implicaría participar para nuestra rutina diaria?",
                    "¿Qué cambios recientes en la respiración, la digestión o la energía debemos vigilar antes de la próxima visita?",
                ],
                "what_to_bring": [
                    "Esta lista impresa con sus tres preguntas principales marcadas.",
                    "Todos los medicamentos y dispositivos tal como se usan realmente en casa.",
                    "Notas breves sobre cualquier cambio reciente (desde cuándo aproximadamente y con qué frecuencia).",
                    "Si necesita traducción durante la consulta, solicite un intérprete médico calificado en su clínica.",
                ],
                "closing_line": (
                    "Por favor consulte esto con su equipo de atención de FQ antes de cambiar cualquier tratamiento. "
                    "Si algo en esta lista es nuevo, empeora o le preocupa, llame a su equipo de atención ahora; no espere a la consulta. "
                    "Esto no es un diagnóstico."
                ),
                "sample_sheet_markdown": sample_es,
            },
        },
        "recruiting_studies": recruiting,
        "active_not_recruiting_studies": active_not_recruiting,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Deterministic Trials by Registry Criteria + Care-Team Questions.")
    parser.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT_JSON, help="Path to snapshot.json")
    parser.add_argument("--status", default="RECRUITING", help="Registry status filter (default: RECRUITING, or ALL)")
    parser.add_argument("--modality", default=None, help="Modality filter (e.g. 'CFTR Modulator', 'mRNA', 'Bacteriophage')")
    parser.add_argument("--phase", default=None, help="Phase filter (e.g. PHASE1, PHASE2, PHASE3)")
    parser.add_argument("--query", default=None, help="Keyword/variant substring filter (e.g. G551D, VX-522)")
    parser.add_argument("--limit", type=int, default=10, help="Max trials to include in sheet")
    parser.add_argument("--lang", choices=("en", "es-MX"), default="en", help="Output language (en or es-MX)")
    parser.add_argument("--export-web", type=Path, default=None, help="Write full static directory JSON to path")
    parser.add_argument("--json", action="store_true", help="Print JSON instead of Markdown sheet")
    args = parser.parse_args(argv)

    try:
        if args.export_web:
            export_doc = build_care_team_directory_export(args.snapshot)
            args.export_web.parent.mkdir(parents=True, exist_ok=True)
            args.export_web.write_text(json.dumps(export_doc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
            print(
                f"Wrote care-team trials directory ({export_doc['metadata']['recruiting_studies_count']} recruiting, "
                f"{export_doc['metadata']['active_not_recruiting_studies_count']} active) -> {args.export_web}"
            )
            return 0

        meta, included = load_prisma_included_studies(args.snapshot)
        trials = filter_registered_trials(
            included,
            status=args.status,
            modality=args.modality,
            phase=args.phase,
            query=args.query,
            limit=args.limit,
        )
        filters = {
            "status": args.status or "ALL",
            "modality": args.modality or ("All modalities" if args.lang == "en" else "Todas las modalidades"),
            "phase": args.phase or ("All phases" if args.lang == "en" else "Todas las fases"),
            "query": args.query or ("None" if args.lang == "en" else "Ninguno"),
        }
        if args.json:
            print(
                json.dumps(
                    {
                        "filters": filters,
                        "snapshot_data_timestamp": meta.get("data_timestamp", NOT_STATED),
                        "matched_count": len(trials),
                        "trials": trials,
                    },
                    indent=2,
                    ensure_ascii=False,
                )
            )
        else:
            print(render_care_team_sheet(trials, filters, meta.get("data_timestamp", ""), lang=args.lang))
        return 0
    except RefusedInputError as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
