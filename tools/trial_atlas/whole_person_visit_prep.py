#!/usr/bin/env python3
"""Bilingual (en & es-MX) Whole-Person Care & Equity CF Clinic-Visit Prep Guide.

Grounded in verified CFF 2025-2026 registry priorities (`cffpr-highlights-2025`),
openFDA modulator Structured Product Labels (`NDA 212273`, `NDA 218730`, `NDA 203188`,
`NDA 206038`, `NDA 210491`), ClinicalTrials.gov registry records, and CFF position/
guideline publications across five clinical domains:
  1. hepatobiliary_dili        - Hepatobiliary surveillance and modulator liver safety
  2. newborn_screening_equity  - Newborn screening timeliness, sweat test access, equity
  3. non_sputum_infection_ntm  - Airway infection, NTM (12%), and non-sputum diagnostics
  4. gi_cfrd_nutrition         - GI/colonoscopy screening, CFRD, and modulator-era BMI
  5. mental_pain_repro_telehealth - Mental health, pain, reproductive health, care team

Safety & Grounding Rules (`cf-ai-safe-use`, `cf-visit-prep`, `AH-4`):
- 100% pure Python 3.10+ standard library, zero network calls, zero LLM calls.
- Refuses personal health information (PHI) via `care_team_trials.validate_query_inputs`.
- Zero medication doses, schedules, or personal clinical eligibility claims.
- Every English output passes `cf-ai-safe-use` (`check_safe_output.py`), with `[T1]`/`[T2]`
  epistemic tier tags, openFDA SPL Boxed Warning isolation (`AH-4`), NLM disclaimer,
  and mandatory urgent-symptom (`(911 in the US)`) and care-team closing lines.

Usage:
    python tools/trial_atlas/whole_person_visit_prep.py --domain all --lang en
    python tools/trial_atlas/whole_person_visit_prep.py --domain hepatobiliary_dili --lang es-MX
    python tools/trial_atlas/whole_person_visit_prep.py --check
    python tools/trial_atlas/whole_person_visit_prep.py --falsify
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from care_team_trials import RefusedInputError, validate_query_inputs  # noqa: E402

NLM_DISCLAIMER = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)
URGENT_LINE_EN = (
    "If a symptom is new, severe, getting worse fast or worrying, above all in a baby or child: "
    "Do not wait for this sheet. Call your CF team's urgent or out-of-hours line, or local "
    "emergency services (911 in the US), now."
)
URGENT_LINE_ES = (
    "Si presenta un síntoma nuevo, grave, que empeora rápido o que le preocupa —especialmente en un bebé o niño—, "
    "no espere a la consulta: llame ahora a la línea de urgencias de su equipo de FQ o a los servicios "
    "locales de emergencia (911 in the US)."
)
CLOSING_LINE_EN = (
    "Please check this with your CF care team before changing anything. "
    "If anything on this list is new, getting worse, or worries you, call your care team now. "
    "Do not wait for the visit. This is not a diagnosis."
)
CLOSING_LINE_ES = (
    "Por favor consulte esto con su equipo de atención de FQ antes de cambiar cualquier tratamiento. "
    "Si algo en esta lista es nuevo, empeora o le preocupa, llame a su equipo de atención ahora; "
    "no espere a la consulta. Esto no es un diagnóstico."
)

FORBIDDEN_UNSAFE_PATTERNS: tuple[tuple[str, re.Pattern[str]], ...] = (
    (
        "clinical_directive",
        re.compile(r"\b(?:you\s+should\s+take|stop\s+taking\s+your|you\s+are\s+eligible|you\s+qualify)\b", re.I),
    ),
    (
        "ungrounded_boxed_extrapolation",
        re.compile(r"\b(?:kalydeco|orkambi|symdeko)\b[^.\n]*?\bboxed\s+warning\b(?![^.\n]*\b(?:not\s+stated|none)\b)", re.I),
    ),
    (
        "clinical_overclaim",
        re.compile(r"\b(?:proven\s+to\s+cure|guarantees\s+clinical\s+benefit|cures\s+cystic\s+fibrosis)\b", re.I),
    ),
)


class VisitPrepVerificationError(ValueError):
    """Raised when safety, provenance, or SPL Boxed Warning isolation fails."""


DOMAIN_ORDER: tuple[str, ...] = (
    "hepatobiliary_dili",
    "newborn_screening_equity",
    "non_sputum_infection_ntm",
    "gi_cfrd_nutrition",
    "mental_pain_repro_telehealth",
)

DOMAIN_CONTENT: dict[str, dict[str, Any]] = {
    "hepatobiliary_dili": {
        "citations": (
            "cff-research-priorities",
            "debray-2026-hepatolcommun",
            "openfda-label",
            "DOI 10.1097/HC9.0000000000001061 (Debray et al. 2026)",
            "openFDA SPL NDA 212273 (TRIKAFTA), NDA 218730 (ALYFTREK), NDA 203188 (KALYDECO), NDA 206038 (ORKAMBI), NDA 210491 (SYMDEKO)",
        ),
        "en": {
            "title": "1. Hepatobiliary Surveillance & Modulator Liver Monitoring [T1]",
            "summary": (
                "[T1] OpenFDA Structured Product Labels report a Boxed Warning for drug-induced liver injury "
                "and liver failure on TRIKAFTA (NDA212273) and ALYFTREK (NDA218730), noting baseline and periodic "
                "hepatic lab monitoring (ALT, AST, and bilirubin). For KALYDECO (NDA203188), ORKAMBI (NDA206038), "
                "and SYMDEKO (NDA210491), a Boxed Warning is not stated in the SPL, while hepatic lab checks remain "
                "part of prescribing labeling. Debray et al. (2026, DOI 10.1097/HC9.0000000000001061) outline CFF "
                "hepatobiliary research gaps including imaging and biomarker tools."
            ),
            "questions": (
                "When are our next routine liver function blood checks (ALT, AST, bilirubin) or liver ultrasound / elastography due?",
                "Are any new over-the-counter supplements, herbal products, or medicines on our home list known to stress the liver?",
                "What early signs (such as dark urine, yellowing eyes, upper-right belly discomfort, or new fatigue) mean we should call the clinic right away?",
            ),
        },
        "es-MX": {
            "title": "1. Vigilancia Hepatobiliar y Monitoreo Hepático con Moduladores [T1]",
            "summary": (
                "[T1] Las fichas oficiales openFDA SPL reportan una Advertencia de Recuadro sobre lesión hepática e "
                "insuficiencia hepática para TRIKAFTA (NDA212273) y ALYFTREK (NDA218730), con monitoreo basal y periódico "
                "de enzimas hepáticas (ALT, AST y bilirrubina). Para KALYDECO (NDA203188), ORKAMBI (NDA206038) y SYMDEKO "
                "(NDA210491), la Advertencia de Recuadro figura como no declarada (not stated) en la etiqueta SPL, aunque se "
                "indican análisis periódicos de transaminasas. Debray et al. (2026, DOI 10.1097/HC9.0000000000001061) "
                "resumen las prioridades de investigación hepatobiliar de la CFF."
            ),
            "questions": (
                "¿Cuándo corresponden los próximos análisis hepáticos (ALT, AST y bilirrubina) o estudios de imagen/elastografía?",
                "¿Alguno de los suplementos o medicamentos adicionales en casa puede sobrecargar el hígado?",
                "¿Qué señales tempranas (orina oscura, ojos amarillos o dolor abdominal derecho) requieren llamar de inmediato a la clínica?",
            ),
        },
    },
    "newborn_screening_equity": {
        "citations": (
            "cff-research-priorities",
            "cffpr-highlights-2025",
            "mcgarry-2025-ijns",
            "DOI 10.3390/ijns11020024 (McGarry et al. 2025)",
            "DOI 10.1001/jamapediatrics.2026.4593 (Ren, Sontag & Faro 2026)",
            "DOI 10.1016/S2213-2600(26)00159-1 (Zemanick et al. 2026)",
            "DOI 10.1016/j.chest.2026.09.107 (Milinic et al. 2026)",
        ),
        "en": {
            "title": "2. Newborn Screening Timeliness, Sweat Testing & Community Equity [T1]",
            "summary": (
                "[T1] The 2025 CFF Patient Registry Highlights document rising diversity among newly diagnosed "
                "individuals: 13.0% Hispanic (versus 10.2% of the total Registry), 3.7% Black, 3.6% multiracial, and "
                "2.1% Asian. CFF newborn screening guidance (McGarry et al. 2025, DOI 10.3390/ijns11020024; Ren, Sontag & "
                "Faro 2026, DOI 10.1001/jamapediatrics.2026.4593) emphasizes pairing IRT panels with broad CFTR variant "
                "sequencing and prompt sweat chloride testing across all ancestral backgrounds, while Zemanick et al. "
                "(2026, DOI 10.1016/S2213-2600(26)00159-1) review sweat chloride interpretation and Milinic et al. "
                "(2026, DOI 10.1016/j.chest.2026.09.107) report lived equity barriers."
            ),
            "questions": (
                "If our baby had an elevated IRT newborn screen or an inconclusive result, how quickly can diagnostic sweat chloride testing and full CFTR sequencing happen?",
                "Can our clinic provide a qualified medical interpreter in our preferred language for genetic counseling and care visits?",
                "How should we interpret intermediate sweat chloride results or follow CRMS/CFSPID guidance over time?",
            ),
        },
        "es-MX": {
            "title": "2. Equidad en Tamiz Neonatal, Prueba de Sudor y Secuenciación Genética [T1]",
            "summary": (
                "[T1] El reporte CFF Patient Registry Highlights 2025 registra mayor diversidad en nuevos diagnósticos: "
                "13.0% hispanos (frente al 10.2% del Registro total), 3.7% negros, 3.6% multirraciales y 2.1% asiáticos. "
                "Las guías de tamiz neonatal de la CFF (McGarry et al. 2025, DOI 10.3390/ijns11020024; Ren, Sontag & Faro "
                "2026, DOI 10.1001/jamapediatrics.2026.4593) recomiendan ampliar la secuenciación genética de CFTR y el acceso "
                "oportuno a la prueba de cloruro en sudor en todas las poblaciones, junto con Zemanick et al. (2026, "
                "DOI 10.1016/S2213-2600(26)00159-1) y Milinic et al. (2026, DOI 10.1016/j.chest.2026.09.107)."
            ),
            "questions": (
                "Si el tamiz neonatal mostró un resultado elevado o no concluyente, ¿qué tan pronto podemos completar la prueba de cloruro en sudor y la secuenciación completa de CFTR?",
                "¿Podemos solicitar un intérprete médico calificado en nuestro idioma preferido para revisar los resultados genéticos?",
                "¿Cómo interpreta el equipo de FQ un resultado intermedio de sudor o el seguimiento tipo CRMS/CFSPID?",
            ),
        },
    },
    "non_sputum_infection_ntm": {
        "citations": (
            "cff-research-priorities",
            "cffpr-highlights-2025",
            "NCT07312734 (SEND-CF observational diagnostics)",
            "NCT06998043 (BX004)",
            "NCT07275905 (AchromoPhage)",
            "NCT04294043 (ABATE)",
            "DOI 10.1099/jmm.0.002150",
        ),
        "en": {
            "title": "3. Airway Infection Surveillance, NTM (12%) & Non-Sputum Diagnostics [T1]",
            "summary": (
                "[T1] In the 2025 CFF Patient Registry Highlights, Pseudomonas aeruginosa prevalence was 22%, MRSA was "
                "12%, and nontuberculous mycobacteria (NTM) was 12% (median first detected age 26 years). Because daily "
                "sputum production decreases for many individuals taking CFTR modulators, observational study NCT07312734 "
                "(SEND-CF) evaluates non-sputum pathogen detection methods alongside sputum culture, while registered phage "
                "and gallium studies (NCT06998043 BX004, NCT07275905 AchromoPhage, NCT04294043 ABATE) remain active "
                "research tracks."
            ),
            "questions": (
                "If producing expectorated sputum has become difficult at clinic visits, which alternative sample types (such as induced sputum, swabs, or diagnostic protocols) does our center use for Pseudomonas and NTM checks?",
                "How frequently does our care team screen respiratory samples for NTM, MRSA, or Achromobacter?",
                "If chronic bacterial colonization persists, how does the team weigh standard inhaled antibiotics versus clinical research options listed in the registry?",
            ),
        },
        "es-MX": {
            "title": "3. Vigilancia de Infecciones, MNT (12%) y Diagnósticos Sin Esputo [T1]",
            "summary": (
                "[T1] En el reporte CFF Patient Registry Highlights 2025, la prevalencia de Pseudomonas aeruginosa fue de "
                "22%, SARM (MRSA) de 12% y micobacterias no tuberculosas (MNT/NTM) de 12% (edad mediana de primera detección "
                "de 26 años). Dado que los moduladores reducen la producción de esputo, el estudio observacional NCT07312734 "
                "(SEND-CF) evalúa métodos diagnósticos alternativos, mientras avanzan estudios registrados con bacteriófagos y "
                "galio intravenoso (NCT06998043 BX004, NCT07275905 AchromoPhage, NCT04294043 ABATE)."
            ),
            "questions": (
                "Si ya casi no producimos esputo en consulta, ¿qué muestras alternativas (esputo inducido, hisopado o protocolos diagnósticos) recomienda el centro para vigilar Pseudomonas y MNT?",
                "¿Con qué frecuencia solicita el equipo cultivos de vigilancia para MNT, SARM o Achromobacter?",
                "Si persiste una infección crónica, ¿qué opciones estándar o ensayos clínicos registrados vale la pena revisar juntos?",
            ),
        },
    },
    "gi_cfrd_nutrition": {
        "citations": ("cff-research-priorities", "cffpr-highlights-2025"),
        "en": {
            "title": "4. Colorectal Cancer Screening, CFRD (24%) & Modulator-Era Nutrition [T1]",
            "summary": (
                "[T1] The 2025 CFF Patient Registry Highlights note that 48% of eligible individuals underwent colonoscopy "
                "screening between 2021 and 2025 (3,891 screened, with 66 colon cancer cases found; screening is noted "
                "starting at age 40 or age 30 after transplant). Additionally, 24% of individuals aged 10 or older had "
                "CF-related diabetes (plus 2% recorded with type 1 or type 2 diabetes), while 29% of adults in 2025 were "
                "overweight and 14% obese (43% combined)."
            ),
            "questions": (
                "When is our next oral glucose tolerance test (OGTT) or CFRD screening, and how are continuous glucose patterns interpreted?",
                "Based on age and transplant history, are colorectal cancer screening (colonoscopy) dates and bowel-prep instructions up to date?",
                "How should nutrition goals and pancreatic enzyme or lipid-soluble vitamin discussions adapt alongside modulator therapy?",
            ),
        },
        "es-MX": {
            "title": "4. Tamizaje Colorrectal, Diabetes Relacionada con FQ (24%) y Nutrición [T1]",
            "summary": (
                "[T1] El reporte CFF Patient Registry Highlights 2025 señala una tasa de colonoscopia preventiva del 48% "
                "entre 2021 y 2025 (3,891 personas evaluadas, con 66 casos de cáncer de colon detectados; el tamizaje inicia a "
                "los 40 años o a los 30 años tras un trasplante). Además, el 24% de las personas de 10 años o más presentaron "
                "diabetes relacionada con FQ (y un 2% adicional diabetes tipo 1 o 2), mientras que en 2025 el 29% de los "
                "adultos tenía sobrepeso y el 14% obesidad (43% combinado)."
            ),
            "questions": (
                "¿Cuándo corresponde la próxima prueba de tolerancia a la glucosa (PTGO/OGTT) o revisión de diabetes relacionada con FQ?",
                "Según edad e historial de trasplante, ¿están al día el tamizaje mediante colonoscopia y su preparación intestinal adaptada para FQ?",
                "¿Cómo conviene ajustar las metas nutricionales, enzimas y vitaminas con nuestro nutriólogo de FQ?",
            ),
        },
    },
    "mental_pain_repro_telehealth": {
        "citations": (
            "cff-research-priorities",
            "cffpr-highlights-2025",
            "DOI 10.1016/j.jcf.2024.11.006 (Bathgate et al. 2025)",
            "DOI 10.1016/j.jcf.2024.11.012 (Dellon et al. 2025)",
            "DOI 10.1016/S2213-2600(26)00120-7 (Jain et al. 2026)",
            "DOI 10.1016/j.jcf.2024.08.007 (Goetz et al. 2024)",
            "DOI 10.1016/j.jcf.2024.09.011 (Brown et al. 2024)",
        ),
        "en": {
            "title": "5. Mental Health, Pain, Reproductive Planning & Team Care Model [T1]",
            "summary": (
                "[T1] The 2025 CFF Patient Registry Highlights record reported depression in 26% and anxiety in 29% of "
                "individuals aged 12 or older, alongside 680 pregnancies, 3,337 people aged 50 or older, and 11,261 "
                "telehealth visits. Recent CFF consensus papers address mental health with modulator use (Bathgate et al. "
                "2025, DOI 10.1016/j.jcf.2024.11.006), pain management (Dellon et al. 2025, DOI 10.1016/j.jcf.2024.11.012), "
                "reproductive health (Jain et al. 2026, DOI 10.1016/S2213-2600(26)00120-7), and evolving telehealth and "
                "care-team coordination (Goetz et al. 2024, DOI 10.1016/j.jcf.2024.08.007; Brown et al. 2024, "
                "DOI 10.1016/j.jcf.2024.09.011)."
            ),
            "questions": (
                "Which sleep, mood, anxiety, or chronic pain symptoms would be helpful to review with the clinic mental health coordinator or social worker?",
                "For family planning, contraception, pregnancy, or fertility questions, can our CF center coordinate pre-conception counseling with obstetrics/genetics specialists?",
                "How can telehealth check-ins, home spirometry, and specialist referrals fit our family's work, school, and travel schedule?",
            ),
        },
        "es-MX": {
            "title": "5. Salud Mental, Dolor, Salud Reproductiva y Modelo de Atención [T1]",
            "summary": (
                "[T1] En el reporte CFF Patient Registry Highlights 2025, entre personas de 12 años o más se registró 26% de "
                "depresión y 29% de ansiedad, además de 680 embarazos, 3,337 personas de 50 años o más y 11,261 consultas por "
                "telemedicina. Documentos de consenso recientes de la CFF abordan salud mental (Bathgate et al. 2025, "
                "DOI 10.1016/j.jcf.2024.11.006), manejo del dolor (Dellon et al. 2025, DOI 10.1016/j.jcf.2024.11.012), "
                "salud reproductiva (Jain et al. 2026, DOI 10.1016/S2213-2600(26)00120-7) y el modelo multidisciplinario con "
                "telemedicina (Goetz et al. 2024, DOI 10.1016/j.jcf.2024.08.007; Brown et al. 2024, "
                "DOI 10.1016/j.jcf.2024.09.011)."
            ),
            "questions": (
                "¿Qué cambios en el sueño, estado de ánimo, ansiedad o dolor crónico conviene revisar con el coordinador de salud mental o trabajo social?",
                "Si tenemos dudas sobre planificación familiar, anticoncepción o embarazo, ¿cómo coordina la clínica la atención conjunta con obstetricia y genética?",
                "¿Qué controles pueden realizarse por telemedicina o espirometría en casa según nuestra estabilidad clínica?",
            ),
        },
    },
}


def _iter_packet_texts(packet: dict[str, Any]) -> list[str]:
    """Collect every user-facing text field individually for strict per-string safety audits."""
    texts: list[str] = [
        str(packet.get("clinical_boundary") or ""),
        str(packet.get("urgent_symptom_notice") or ""),
        str(packet.get("closing_notice") or ""),
        str(packet.get("nlm_disclaimer") or ""),
    ]
    for sec in packet.get("sections", []):
        texts.append(str(sec.get("title") or ""))
        texts.append(str(sec.get("summary") or ""))
        for q in sec.get("care_team_questions", []):
            texts.append(str(q))
    return texts


def verify_visit_prep_packet(packet: dict[str, Any]) -> None:
    """Verify non-clinical boundaries, urgent notice, NLM disclaimer, and AH-4 Boxed Warning rules."""
    for text_piece in _iter_packet_texts(packet):
        for rule_name, pat in FORBIDDEN_UNSAFE_PATTERNS:
            if pat.search(text_piece):
                raise VisitPrepVerificationError(f"Safety rule failed ({rule_name}) in visit-prep output.")
    serialized = json.dumps(packet, ensure_ascii=False)
    for rule_name, pat in FORBIDDEN_UNSAFE_PATTERNS:
        if pat.search(serialized):
            raise VisitPrepVerificationError(f"Safety rule failed ({rule_name}) in visit-prep output.")
    if "(911 in the US)" not in str(packet.get("urgent_symptom_notice") or ""):
        raise VisitPrepVerificationError("Missing mandatory '(911 in the US)' urgent symptom notice.")
    if "does not evaluate safety or efficacy" not in str(packet.get("nlm_disclaimer") or ""):
        raise VisitPrepVerificationError("Missing mandatory NLM disclaimer.")
    boundary = str(packet.get("clinical_boundary") or "").lower()
    if "not" not in boundary and "no " not in boundary:
        raise VisitPrepVerificationError("Missing explicit non-clinical boundary statement.")


def build_whole_person_packet(
    domain: str = "all",
    lang: str = "en",
    custom_note_query: str | None = None,
) -> dict[str, Any]:
    """Build a verified bilingual whole-person clinic-visit preparation packet."""
    validate_query_inputs(domain, lang, custom_note_query)
    dom_norm = (domain or "all").strip().lower()
    if dom_norm != "all" and dom_norm not in DOMAIN_CONTENT:
        raise ValueError(f"Unknown --domain {domain!r}. Allowed: 'all' or {list(DOMAIN_ORDER)}")

    is_es = (lang or "en").strip().lower().startswith("es")
    lang_key = "es-MX" if is_es else "en"
    selected = list(DOMAIN_ORDER) if dom_norm == "all" else [dom_norm]

    boundary = (
        "Solo para preparación informativa de consulta (Tier T1). No es consejo médico ni diagnóstico; "
        "verifique los códigos NCT en el sitio web oficial con su equipo de atención de FQ."
        if is_es
        else (
            "Educational clinic-visit preparation only (Tier T1). Independent open-source research, "
            "not medical advice or diagnosis. Verify NCT records on the official website with your CF care team."
        )
    )

    sections: list[dict[str, Any]] = []
    for key in selected:
        entry = DOMAIN_CONTENT[key]
        localized = entry[lang_key]
        sections.append(
            {
                "domain_id": key,
                "title": localized["title"],
                "summary": localized["summary"],
                "care_team_questions": list(localized["questions"]),
                "primary_citations": list(entry["citations"]),
            }
        )

    packet: dict[str, Any] = {
        "schema_version": 1,
        "tool": "whole_person_visit_prep",
        "lang": lang_key,
        "epistemic_tier": "T1",
        "source_catalog_ids": [
            "cff-research-priorities",
            "cffpr-highlights-2025",
            "debray-2026-hepatolcommun",
            "mcgarry-2025-ijns",
            "openfda-label",
            "clinicaltrials-gov",
        ],
        "clinical_boundary": boundary,
        "urgent_symptom_notice": URGENT_LINE_ES if is_es else URGENT_LINE_EN,
        "closing_notice": CLOSING_LINE_ES if is_es else CLOSING_LINE_EN,
        "nlm_disclaimer": NLM_DISCLAIMER,
        "sections": sections,
    }
    verify_visit_prep_packet(packet)
    return packet


def render_whole_person_markdown(packet: dict[str, Any]) -> str:
    """Render a whole-person packet into Markdown suitable for `check_safe_output.py` auditing."""
    is_es = packet.get("lang") == "es-MX"
    heading = (
        "# Guía Bilingüe de Preparación para Consulta de FQ (Prioridades CFF 2025-2026) [T1]"
        if is_es
        else "# CF Whole-Person Clinic Visit Prep Guide (CFF 2025-2026 Priorities) [T1]"
    )
    lines = [
        heading,
        "",
        f"{packet['clinical_boundary']} {packet['urgent_symptom_notice']}",
        "",
        f"[T1] Please check them on the official website with your CF care team. {packet['nlm_disclaimer']}",
        "",
    ]
    for sec in packet["sections"]:
        lines.append(f"## {sec['title']}")
        lines.append(sec["summary"])
        lines.append("")
        lines.append("**Preguntas sugeridas para su equipo de FQ:**" if is_es else "**Suggested Questions for Your CF Care Team:**")
        for q in sec["care_team_questions"]:
            lines.append(f"- {q}")
        lines.append("")

    lines.append(packet["closing_notice"])
    return "\n".join(lines)


def run_falsification_suite() -> dict[str, Any]:
    """Run deterministic negative controls (--falsify) covering PHI, AH-4 Boxed Warning, and directives."""
    controls: list[dict[str, Any]] = []

    # Control 1: PHI query refused
    try:
        build_whole_person_packet(domain="all", custom_note_query="DOB 02/14/2018 MRN 55119988")
        controls.append({"control": "phi_input_refusal", "passed": False, "detail": "PHI input was accepted"})
    except RefusedInputError as exc:
        controls.append({"control": "phi_input_refusal", "passed": True, "detail": str(exc)})

    # Control 2: Extrapolated Boxed Warning on Kalydeco without 'not stated' fails AH-4 guard even across multi-domain all-packet serializations
    poisoned_bw = build_whole_person_packet(domain="all", lang="en")
    poisoned_bw["sections"][0]["summary"] = "KALYDECO carries a Boxed Warning requiring immediate dose discontinuation."
    poisoned_bw["sections"][1]["summary"] = "Restrictions: none."
    try:
        verify_visit_prep_packet(poisoned_bw)
        controls.append({"control": "extrapolated_boxed_warning_blocked", "passed": False, "detail": "AH-4 violation slipped past"})
    except VisitPrepVerificationError as exc:
        controls.append({"control": "extrapolated_boxed_warning_blocked", "passed": True, "detail": str(exc)})

    # Control 3: Directive eligibility claim fails safety guard
    poisoned_dir = build_whole_person_packet(domain="non_sputum_infection_ntm", lang="en")
    poisoned_dir["sections"][0]["summary"] = "You qualify immediately and should take this antibiotic."
    try:
        verify_visit_prep_packet(poisoned_dir)
        controls.append({"control": "directive_claim_blocked", "passed": False, "detail": "Directive claim slipped past"})
    except VisitPrepVerificationError as exc:
        controls.append({"control": "directive_claim_blocked", "passed": True, "detail": str(exc)})

    all_passed = all(c["passed"] for c in controls)
    return {
        "tool": "whole_person_visit_prep",
        "falsify_passed": all_passed,
        "negative_controls": controls,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--domain",
        type=str,
        choices=("all", *DOMAIN_ORDER),
        default="all",
        help=f"Topic domain: all or {','.join(DOMAIN_ORDER)}",
    )
    parser.add_argument("--lang", type=str, choices=("en", "es-MX"), default="en", help="Output language")
    parser.add_argument("--json", action="store_true", help="Emit JSON packet")
    parser.add_argument("--check", action="store_true", help="Run fast deterministic self-check")
    parser.add_argument("--falsify", action="store_true", help="Run negative-control falsification suite")
    args = parser.parse_args(argv)

    if args.falsify:
        f_rep = run_falsification_suite()
        if args.json:
            print(json.dumps(f_rep, indent=2))
        else:
            print(f"PASS falsify whole_person_visit_prep ({len(f_rep['negative_controls'])}/3 negative controls)" if f_rep["falsify_passed"] else "FAIL falsify whole_person_visit_prep")
        return 0 if f_rep["falsify_passed"] else 1

    if args.check:
        pkt_en = build_whole_person_packet(domain="all", lang="en")
        pkt_es = build_whole_person_packet(domain="all", lang="es-MX")
        md_en = render_whole_person_markdown(pkt_en)
        f_rep = run_falsification_suite()
        if not f_rep["falsify_passed"] or len(pkt_en["sections"]) != 5 or len(pkt_es["sections"]) != 5 or "(911 in the US)" not in md_en:
            print("FAIL whole_person_visit_prep self-check", file=sys.stderr)
            return 1
        print(
            f"PASS whole_person_visit_prep: domains={len(pkt_en['sections'])} "
            f"langs=en,es-MX falsify_controls={len(f_rep['negative_controls'])}/3"
        )
        return 0

    try:
        pkt = build_whole_person_packet(domain=args.domain, lang=args.lang)
    except (RefusedInputError, VisitPrepVerificationError, ValueError) as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(pkt, indent=2, ensure_ascii=False))
    else:
        print(render_whole_person_markdown(pkt))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
