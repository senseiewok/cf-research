#!/usr/bin/env python3
"""Deterministic 2026 CFF Priority Clinical Research Radar (T1 Registry Watchlist).

Tracks Cystic Fibrosis Foundation (CFF) 2025-2026 priority pipeline studies across five
pillars: genetic therapies, next-generation modulators, observational trial readiness,
infection research, and modulator-era care-model optimization.

Strict Invariants (`cf-research-context`, `security-browsing`, `AH-1`..`AH-4`):
- 100% pure Python 3.10+ standard library, zero network calls in offline checks/tests.
- Strictly partitions PRISMA-included interventional trials (`940`) and raw snapshot
  interventional studies (`1,053` in `tools/trial_atlas/snapshot.json`) from
  `OBSERVATIONAL` / post-snapshot watchlist studies (`NCT06504589` REACH and
  `NCT07312734` SEND-CF) so interventional denominators are never contaminated.
- Quarantines any external registry text inside `<untrusted_data_NONCE>` cages with
  Trojan-Source bidi escaping and homoglyph breakout neutralization.
- Enforces the verbatim National Library of Medicine (NLM) disclaimer and non-clinical
  research boundary without making personal eligibility promises.

Usage:
    python tools/trial_atlas/cff_priority_radar.py --check
    python tools/trial_atlas/cff_priority_radar.py --falsify
    python tools/trial_atlas/cff_priority_radar.py --json
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

from care_team_trials import (  # noqa: E402
    DEFAULT_SNAPSHOT_JSON,
    NOT_STATED,
    load_prisma_included_studies,
)

NLM_DISCLAIMER = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)
CLINICAL_BOUNDARY = (
    "Independent open-source research overview only (Tier T1 registry metadata); "
    "not medical advice and never an assessment of personal trial eligibility."
)
ALLOWED_REGISTRY_HOSTS: frozenset[str] = frozenset(
    {"clinicaltrials.gov", "www.clinicaltrials.gov"}
)

TROJAN_BIDI_CODEPOINTS: frozenset[int] = frozenset(
    {
        0x200B, 0x200C, 0x200D, 0x200E, 0x200F,
        0x202A, 0x202B, 0x202C, 0x202D, 0x202E,
        0x2060, 0x2066, 0x2067, 0x2068, 0x2069,
        0xFEFF,
    }
)
CONFUSABLE_ASCII_MAP = str.maketrans(
    {
        "\u0430": "a", "\u0435": "e", "\u043e": "o", "\u0440": "p",
        "\u0441": "c", "\u0443": "y", "\u0445": "x", "\u0501": "d",
        "\u0456": "i", "\u0455": "s", "\u0442": "t", "\u043d": "n",
        "\u0410": "A", "\u0412": "B", "\u0421": "C", "\u0500": "D",
        "\u0415": "E", "\u041d": "H", "\u0406": "I", "\u041a": "K",
        "\u041c": "M", "\u041e": "O", "\u0420": "P", "\u0405": "S",
        "\u0422": "T", "\u0425": "X", "\u0423": "Y",
    }
)
SOURCE_CATALOG_IDS: tuple[str, ...] = (
    "cff-research-priorities",
    "cffpr-highlights-2025",
    "clinicaltrials-gov",
)
UNTRUSTED_CLOSE_LOOKALIKE_RE = re.compile(r"<\s*/\s*[^\s>]*data[^\s>]*>", re.I)

FORBIDDEN_ELIGIBILITY_PHRASES: tuple[re.Pattern[str], ...] = (
    re.compile(r"\b(?:you\s+are\s+eligible|you\s+may\s+be\s+eligible|you\s+qualify|qualifies\s+for)\b", re.I),
    re.compile(r"\b(?:proven\s+to\s+cure|guarantees\s+clinical\s+benefit|safe\s+and\s+effective\s+for\s+patients)\b", re.I),
)


class PriorityRadarError(ValueError):
    """Raised when watchlist cohort separation, provenance, or quarantine invariants fail."""


# Verified 18-study CFF 2025-2026 priority watchlist.
# Note: NCT04135495 is ELX-02 (G542X readthrough), whereas NCT04596319 is Armata AP-PA02 phage,
# and NCT06515002/NCT06962852 are BI 3720931 lentiviral trials (NCT05033080 is SKYLINE 102 VTD).
CFF_2026_PRIORITY_WATCHLIST: tuple[dict[str, Any], ...] = (
    {
        "nct_id": "NCT06237335",
        "program_label": "RCT2100 (ReCode inhaled CFTR mRNA-LNP)",
        "priority_pillar": "genetic_therapies",
        "study_type": "INTERVENTIONAL",
        "expected_in_snapshot_2026_10": True,
        "expected_in_prisma_940": True,
        "population_note": "Healthy participants and participants with Cystic Fibrosis (Phase 2).",
        "clinical_readout": "safety, tolerability, and CFTR functional restoration exploration",
    },
    {
        "nct_id": "NCT05248230",
        "program_label": "4D-710 (4DMT AEROW aerosolized A101 AAV CFTRdeltaR)",
        "priority_pillar": "genetic_therapies",
        "study_type": "INTERVENTIONAL",
        "expected_in_snapshot_2026_10": True,
        "expected_in_prisma_940": True,
        "population_note": "Adults with CF ineligible for or not taking CFTR modulators.",
        "clinical_readout": "safety, tolerability, bronchial CFTR expression, and ppFEV1",
    },
    {
        "nct_id": "NCT06429176",
        "program_label": "SPL84 (SpliSense inhaled ASO for 3849+10kb C->T)",
        "priority_pillar": "genetic_therapies",
        "study_type": "INTERVENTIONAL",
        "expected_in_snapshot_2026_10": True,
        "expected_in_prisma_940": True,
        "population_note": "People with CF carrying the 3849+10kb C->T splicing mutation.",
        "clinical_readout": "safety, tolerability, PK, and preliminary ppFEV1 efficacy",
    },
    {
        "nct_id": "NCT06526923",
        "program_label": "SP-101 (Spirovant inhaled AAV2.5T + doxorubicin enhancer)",
        "priority_pillar": "genetic_therapies",
        "study_type": "INTERVENTIONAL",
        "expected_in_snapshot_2026_10": True,
        "expected_in_prisma_940": True,
        "population_note": "Participants with CF (Phase 1/2).",
        "clinical_readout": "safety and tolerability of inhaled AAV2.5T transduction",
    },
    {
        "nct_id": "NCT05504837",
        "program_label": "KB407 (Krystal nebulized HSV-1 full-length CFTR vector)",
        "priority_pillar": "genetic_therapies",
        "study_type": "INTERVENTIONAL",
        "expected_in_snapshot_2026_10": True,
        "expected_in_prisma_940": True,
        "population_note": "Participants with CF independent of CFTR mutation (Phase 1).",
        "clinical_readout": "safety, tolerability, and vector shedding/expression",
    },
    {
        "nct_id": "NCT07485543",
        "program_label": "BMD003 (Inhaled CFTR mRNA candidate)",
        "priority_pillar": "genetic_therapies",
        "study_type": "INTERVENTIONAL",
        "expected_in_snapshot_2026_10": True,
        "expected_in_prisma_940": True,
        "population_note": "Participants aged >= 12 years with Cystic Fibrosis (Early Phase 1).",
        "clinical_readout": "safety and preliminary clinical/functional endpoints",
    },
    {
        "nct_id": "NCT05668741",
        "program_label": "VX-522 (Vertex inhaled CFTR mRNA-LNP historical benchmark)",
        "priority_pillar": "genetic_therapies",
        "study_type": "INTERVENTIONAL",
        "expected_in_snapshot_2026_10": True,
        "expected_in_prisma_940": True,
        "population_note": "Participants aged >= 18 years ineligible for CFTR modulators.",
        "clinical_readout": "safety, tolerability, and delivery lessons for inhaled LNP-mRNA",
    },
    {
        "nct_id": "NCT06515002",
        "program_label": "BI 3720931 (Boehringer Ingelheim Lenticlair 1 lentiviral CFTR)",
        "priority_pillar": "genetic_therapies",
        "study_type": "INTERVENTIONAL",
        "expected_in_snapshot_2026_10": True,
        "expected_in_prisma_940": True,
        "population_note": "Adults with CF not eligible for CFTR modulators (registry status: TERMINATED).",
        "clinical_readout": "tolerability, safety, and ppFEV1 change",
    },
    {
        "nct_id": "NCT06962852",
        "program_label": "BI 3720931 (Boehringer Ingelheim Lenticlair-ON safety follow-up)",
        "priority_pillar": "genetic_therapies",
        "study_type": "INTERVENTIONAL",
        "expected_in_snapshot_2026_10": True,
        "expected_in_prisma_940": True,
        "population_note": "Long-term follow-up of participants previously dosed with BI 3720931.",
        "clinical_readout": "long-term safety monitoring following lentiviral vector dosing",
    },
    {
        "nct_id": "NCT04135495",
        "program_label": "ELX-02 (Class I G542X nonsense readthrough benchmark)",
        "priority_pillar": "genetic_therapies",
        "study_type": "INTERVENTIONAL",
        "expected_in_snapshot_2026_10": True,
        "expected_in_prisma_940": True,
        "population_note": "Participants carrying at least one G542X nonsense allele (Phase 2 COMPLETED).",
        "clinical_readout": "safety, PK/PD, and sweat chloride exploratory readthrough response",
    },
    {
        "nct_id": "NCT07283770",
        "program_label": "VX-581 (Vertex next-generation modulator Phase 1 PK/safety)",
        "priority_pillar": "next_gen_modulators",
        "study_type": "INTERVENTIONAL",
        "expected_in_snapshot_2026_10": True,
        "expected_in_prisma_940": True,
        "population_note": "Healthy adult participants only (Phase 1 PK/safety; does not enroll people with CF).",
        "clinical_readout": "dose-escalation safety, tolerability, and pharmacokinetics with TEZ/D-IVA",
    },
    {
        "nct_id": "NCT06504589",
        "program_label": "REACH (CFF TDN Observational Trial-Readiness Cohort)",
        "priority_pillar": "observational_readiness",
        "study_type": "OBSERVATIONAL",
        "expected_in_snapshot_2026_10": False,
        "expected_in_prisma_940": False,
        "population_note": "People with CF who are ineligible for or not taking CFTR modulators (external OBSERVATIONAL).",
        "clinical_readout": "longitudinal clinical trajectories and biospecimen readiness (ppFEV1 at 12 months)",
        "registered_status_snapshot_verified": "RECRUITING",
        "registered_brief_title": "A Research Study to Advance the CF Therapeutics Pipeline for People Without Modulators",
        "verified_api_date": "2025-10-01",
    },
    {
        "nct_id": "NCT06998043",
        "program_label": "BX004 (BiomX nebulized bacteriophage for P. aeruginosa)",
        "priority_pillar": "infection_initiative",
        "study_type": "INTERVENTIONAL",
        "expected_in_snapshot_2026_10": True,
        "expected_in_prisma_940": True,
        "population_note": "Participants with CF and chronic Pseudomonas aeruginosa pulmonary infection (Phase 2).",
        "clinical_readout": "sputum bacterial density reduction and safety/tolerability",
    },
    {
        "nct_id": "NCT07275905",
        "program_label": "AchromoPhage (Bacteriophage cocktail for chronic Achromobacter)",
        "priority_pillar": "infection_initiative",
        "study_type": "INTERVENTIONAL",
        "expected_in_snapshot_2026_10": True,
        "expected_in_prisma_940": True,
        "population_note": "Adults with CF and chronic Achromobacter airway infection (Phase 1).",
        "clinical_readout": "safety, tolerability, and microbiological kinetics across inhaled/IV routes",
    },
    {
        "nct_id": "NCT04596319",
        "program_label": "AP-PA02 (Armata inhaled P. aeruginosa bacteriophage benchmark)",
        "priority_pillar": "infection_initiative",
        "study_type": "INTERVENTIONAL",
        "expected_in_snapshot_2026_10": True,
        "expected_in_prisma_940": True,
        "population_note": "Participants with CF and chronic P. aeruginosa lung infection (Phase 1/2 COMPLETED).",
        "clinical_readout": "safety, tolerability, and phage/bacterial recovery",
    },
    {
        "nct_id": "NCT04294043",
        "program_label": "ABATE (IV gallium nitrate for NTM pulmonary infection)",
        "priority_pillar": "infection_initiative",
        "study_type": "INTERVENTIONAL",
        "expected_in_snapshot_2026_10": True,
        "expected_in_prisma_940": False,  # Excluded by PRISMA X4 while present in raw 1,053 snapshot
        "population_note": "Participants with CF and nontuberculous mycobacteria (NTM) lung infection (Phase 1).",
        "clinical_readout": "safety, pharmacokinetics, and NTM microbiological response",
    },
    {
        "nct_id": "NCT07312734",
        "program_label": "SEND-CF (Non-sputum pulmonary microbial diagnostics vs. sputum culture)",
        "priority_pillar": "infection_initiative",
        "study_type": "OBSERVATIONAL",
        "expected_in_snapshot_2026_10": False,
        "expected_in_prisma_940": False,
        "population_note": "People with CF evaluated for alternative non-sputum pathogen detection (external OBSERVATIONAL).",
        "clinical_readout": "sensitivity/specificity of alternative microbial diagnostics vs. sputum culture",
        "registered_status_snapshot_verified": "RECRUITING",
        "registered_brief_title": "Study to Enable New Diagnostics for Pulmonary Microbes in People With CF",
        "verified_api_date": "2026-07-02",
    },
    {
        "nct_id": "NCT06654752",
        "program_label": "STOP-PEDS (Streamlined Treatment of Pulmonary Exacerbations in Pediatrics)",
        "priority_pillar": "care_model",
        "study_type": "INTERVENTIONAL",
        "expected_in_snapshot_2026_10": True,
        "expected_in_prisma_940": True,
        "population_note": "Pediatric CF pulmonary exacerbations comparing immediate vs. tailored oral antibiotics.",
        "clinical_readout": "comparative exacerbation recovery and antibiotic stewardship in pediatric care",
    },
)


def is_allowed_registry_url(raw_url: str) -> bool:
    """Validate URL scheme and exact IDNA-normalized host against ALLOWED_REGISTRY_HOSTS."""
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
    """Escape Trojan-Source bidi controls and homoglyph tag breakouts inside <untrusted_data_NONCE>."""
    active_nonce = nonce or secrets.token_hex(8)
    escaped_chars: list[str] = []
    trojan_escaped = 0
    for ch in str(raw_text or ""):
        cp = ord(ch)
        cat = unicodedata.category(ch)
        if cp in TROJAN_BIDI_CODEPOINTS or (cat in ("Cf", "Cc") and ch not in ("\n", "\r", "\t")):
            escaped_chars.append(f"<U+{cp:04X}>")
            trojan_escaped += 1
        else:
            escaped_chars.append(ch)
    sanitized = "".join(escaped_chars)

    skeleton = sanitized.translate(CONFUSABLE_ASCII_MAP)
    lookalikes_escaped = 0
    if re.search(r"<\s*/\s*untrusted", skeleton, re.I) or UNTRUSTED_CLOSE_LOOKALIKE_RE.search(skeleton):
        def _neutralize(match: re.Match[str]) -> str:
            nonlocal lookalikes_escaped
            lookalikes_escaped += 1
            return match.group(0).replace("<", "&lt;").replace(">", "&gt;")

        sanitized = re.sub(r"<\s*/\s*[^>]+>", _neutralize, sanitized)

    caged = f"<untrusted_data_{active_nonce}>\n{sanitized}\n</untrusted_data_{active_nonce}>"
    return {
        "nonce": active_nonce,
        "caged_text": caged,
        "trojan_chars_escaped": trojan_escaped,
        "boundary_lookalikes_escaped": lookalikes_escaped,
    }


def verify_radar_report(report: dict[str, Any]) -> None:
    """Enforce strict separation between interventional snapshot cohorts and observational/external rows."""
    serialized = json.dumps(report, ensure_ascii=False)
    for pat in FORBIDDEN_ELIGIBILITY_PHRASES:
        if pat.search(serialized):
            raise PriorityRadarError("Report contains prohibited clinical overclaim or personal eligibility phrasing.")

    if "does not evaluate safety or efficacy" not in str(report.get("nlm_disclaimer") or ""):
        raise PriorityRadarError("Report missing mandatory NLM ClinicalTrials.gov disclaimer.")

    interventional_rows = report.get("interventional_snapshot_watchlist", [])
    external_rows = report.get("observational_or_external_watchlist", [])

    interventional_ncts = {str(r.get("nct_id") or "") for r in interventional_rows}
    for row in interventional_rows:
        if row.get("study_type") == "OBSERVATIONAL":
            raise PriorityRadarError(
                f"Cohort mixing violation: OBSERVATIONAL study {row.get('nct_id')} placed in interventional_snapshot_watchlist."
            )
    for row in external_rows:
        nid = str(row.get("nct_id") or "")
        if nid in interventional_ncts:
            raise PriorityRadarError(
                f"Cohort mixing violation: external/observational study {nid} duplicated in interventional_snapshot_watchlist."
            )
        if row.get("included_in_prisma_940"):
            raise PriorityRadarError(
                f"Cohort mixing violation: external/observational study {nid} claimed included_in_prisma_940=True."
            )
        reg_url = row.get("registry_url")
        if reg_url and not is_allowed_registry_url(str(reg_url)):
            raise PriorityRadarError(
                f"Host validation failed: external study {nid} carries non-canonical registry URL."
            )


def build_priority_radar_report(
    snapshot_path: Path = DEFAULT_SNAPSHOT_JSON,
    watchlist: tuple[dict[str, Any], ...] = CFF_2026_PRIORITY_WATCHLIST,
) -> dict[str, Any]:
    """Cross-verify the 2026 CFF priority watchlist against snapshot.json (1,053 raw / 940 PRISMA)."""
    raw_data = json.loads(snapshot_path.read_text(encoding="utf-8"))
    raw_studies = raw_data.get("studies", [])
    raw_by_nct = {str(s.get("nct_id")): s for s in raw_studies if s.get("nct_id")}

    meta, prisma_studies = load_prisma_included_studies(snapshot_path)
    prisma_ncts = {str(s.get("nct_id")) for s in prisma_studies if s.get("nct_id")}

    interventional_rows: list[dict[str, Any]] = []
    external_rows: list[dict[str, Any]] = []
    errors: list[str] = []

    for entry in watchlist:
        nct_id = str(entry["nct_id"])
        expected_snap = bool(entry.get("expected_in_snapshot_2026_10"))
        snap_rec = raw_by_nct.get(nct_id)
        in_raw = snap_rec is not None
        in_prisma = nct_id in prisma_ncts

        if expected_snap and not in_raw:
            errors.append(f"Expected snapshot trial {nct_id} missing from {snapshot_path.name}")
        if not expected_snap and in_prisma:
            errors.append(f"External/observational trial {nct_id} leaked into PRISMA 940 cohort")
        if expected_snap and bool(entry.get("expected_in_prisma_940")) != in_prisma:
            errors.append(
                f"PRISMA inclusion mismatch for {nct_id}: expected={entry.get('expected_in_prisma_940')} observed={in_prisma}"
            )

        if entry.get("study_type") == "INTERVENTIONAL" and expected_snap:
            brief_title = snap_rec.get("brief_title", NOT_STATED) if snap_rec else NOT_STATED
            overall_status = snap_rec.get("overall_status", NOT_STATED) if snap_rec else NOT_STATED
            phases = list(snap_rec.get("phases") or [NOT_STATED]) if snap_rec else [NOT_STATED]
            lead_sponsor = snap_rec.get("lead_sponsor", NOT_STATED) if snap_rec else NOT_STATED
            interventional_rows.append(
                {
                    "nct_id": nct_id,
                    "program_label": entry["program_label"],
                    "priority_pillar": entry["priority_pillar"],
                    "study_type": "INTERVENTIONAL",
                    "provenance_tier": "T1",
                    "in_raw_snapshot_1053": in_raw,
                    "included_in_prisma_940": in_prisma,
                    "overall_status": overall_status,
                    "phases": phases,
                    "lead_sponsor": lead_sponsor,
                    "brief_title": brief_title,
                    "population_note": entry["population_note"],
                    "clinical_readout": entry["clinical_readout"],
                    "catalog_source_id": "cff-research-priorities",
                    "registry_url": f"https://clinicaltrials.gov/study/{nct_id}",
                }
            )
        else:
            external_rows.append(
                {
                    "nct_id": nct_id,
                    "program_label": entry["program_label"],
                    "priority_pillar": entry["priority_pillar"],
                    "study_type": entry["study_type"],
                    "provenance_tier": "T1",
                    "in_raw_snapshot_1053": in_raw,
                    "included_in_prisma_940": False,
                    "overall_status": entry.get("registered_status_snapshot_verified", NOT_STATED),
                    "brief_title": entry.get("registered_brief_title", NOT_STATED),
                    "verified_api_date": entry.get("verified_api_date", NOT_STATED),
                    "population_note": entry["population_note"],
                    "clinical_readout": entry["clinical_readout"],
                    "catalog_source_id": "cff-research-priorities",
                    "registry_url": f"https://clinicaltrials.gov/study/{nct_id}",
                }
            )

    if errors:
        raise PriorityRadarError("; ".join(errors))

    report: dict[str, Any] = {
        "schema_version": 1,
        "tool": "cff_priority_radar",
        "epistemic_tier": "T1",
        "source_catalog_ids": list(SOURCE_CATALOG_IDS),
        "snapshot_data_timestamp": meta.get("data_timestamp", "2026-10-14T14:12:11.425555+00:00"),
        "snapshot_raw_interventional_count": len(raw_studies),
        "prisma_included_interventional_count": len(prisma_studies),
        "watchlist_total_count": len(watchlist),
        "interventional_snapshot_count": len(interventional_rows),
        "observational_or_external_count": len(external_rows),
        "clinical_boundary": CLINICAL_BOUNDARY,
        "nlm_disclaimer": NLM_DISCLAIMER,
        "interventional_snapshot_watchlist": interventional_rows,
        "observational_or_external_watchlist": external_rows,
    }
    verify_radar_report(report)
    return report


def run_falsification_suite(snapshot_path: Path = DEFAULT_SNAPSHOT_JSON) -> dict[str, Any]:
    """Run deterministic negative controls (--falsify) verifying every poisoned scenario fails closed."""
    controls: list[dict[str, Any]] = []

    # Control 1: Planted missing NCT expected in snapshot.json must raise PriorityRadarError
    poisoned_watchlist = CFF_2026_PRIORITY_WATCHLIST + (
        {
            "nct_id": "NCT99999999",
            "program_label": "Planted Phantom Trial",
            "priority_pillar": "genetic_therapies",
            "study_type": "INTERVENTIONAL",
            "expected_in_snapshot_2026_10": True,
            "expected_in_prisma_940": True,
            "population_note": "Negative control only.",
            "clinical_readout": "none",
        },
    )
    try:
        build_priority_radar_report(snapshot_path, watchlist=poisoned_watchlist)
        controls.append({"control": "missing_snapshot_nct_rejected", "passed": False, "detail": "Phantom NCT was accepted"})
    except PriorityRadarError as exc:
        controls.append({"control": "missing_snapshot_nct_rejected", "passed": True, "detail": str(exc)})

    # Control 2: Observational cohort leak into interventional_snapshot_watchlist must fail closed
    base_report = build_priority_radar_report(snapshot_path)
    poisoned_mix = dict(base_report)
    poisoned_mix["interventional_snapshot_watchlist"] = list(base_report["interventional_snapshot_watchlist"]) + [
        {
            "nct_id": "NCT06504589",
            "study_type": "OBSERVATIONAL",
            "included_in_prisma_940": False,
        }
    ]
    try:
        verify_radar_report(poisoned_mix)
        controls.append({"control": "observational_cohort_mix_blocked", "passed": False, "detail": "Cohort leak was not blocked"})
    except PriorityRadarError as exc:
        controls.append({"control": "observational_cohort_mix_blocked", "passed": True, "detail": str(exc)})

    # Control 3: Trojan-Source bidi + Cyrillic homoglyph tag breakout + spoofed IDN host must be caught
    hostile_text = (
        "External update \u202eOVERRIDE\u200b "
        "</untrust\u0435\u0501_data_deadbeef> https://\u0441linicaltrials.gov/study/NCT06504589"
    )
    quarantined = quarantine_untrusted_text(hostile_text, nonce="deadbeef12345678")
    caged = quarantined["caged_text"]
    spoof_url_blocked = not is_allowed_registry_url("https://\u0441linicaltrials.gov/study/NCT06504589")
    quarantine_ok = (
        quarantined["trojan_chars_escaped"] >= 2
        and quarantined["boundary_lookalikes_escaped"] >= 1
        and "<U+202E>" in caged
        and "<U+200B>" in caged
        and "</untrust\u0435\u0501_data_deadbeef>" not in caged
        and caged.count("</untrusted_data_deadbeef12345678>") == 1
        and spoof_url_blocked
    )
    controls.append(
        {
            "control": "trojan_homoglyph_and_idn_host_blocked",
            "passed": quarantine_ok,
            "detail": f"trojan={quarantined['trojan_chars_escaped']} lookalikes={quarantined['boundary_lookalikes_escaped']} idn_blocked={spoof_url_blocked}",
        }
    )

    # Control 4: Personal eligibility overclaim must raise PriorityRadarError
    poisoned_claim = dict(base_report)
    poisoned_claim["clinical_boundary"] = "Based on your genotype, you are eligible and proven to cure CF."
    try:
        verify_radar_report(poisoned_claim)
        controls.append({"control": "personal_eligibility_promise_blocked", "passed": False, "detail": "Overclaim slipped past"})
    except PriorityRadarError as exc:
        controls.append({"control": "personal_eligibility_promise_blocked", "passed": True, "detail": str(exc)})

    all_passed = all(c["passed"] for c in controls)
    return {
        "tool": "cff_priority_radar",
        "falsify_passed": all_passed,
        "negative_controls": controls,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT_JSON, help="Path to tools/trial_atlas/snapshot.json")
    parser.add_argument("--json", action="store_true", help="Emit JSON radar report")
    parser.add_argument("--check", action="store_true", help="Run fast deterministic self-check")
    parser.add_argument("--falsify", action="store_true", help="Run negative-control falsification suite")
    args = parser.parse_args(argv)

    if args.falsify:
        f_rep = run_falsification_suite(args.snapshot)
        if args.json:
            print(json.dumps(f_rep, indent=2))
        else:
            status_word = "PASS" if f_rep["falsify_passed"] else "FAIL"
            print(f"{status_word} falsify cff_priority_radar ({sum(1 for c in f_rep['negative_controls'] if c['passed'])}/{len(f_rep['negative_controls'])} negative controls)")
        return 0 if f_rep["falsify_passed"] else 1

    report = build_priority_radar_report(args.snapshot)

    if args.check:
        f_rep = run_falsification_suite(args.snapshot)
        if not f_rep["falsify_passed"] or report["watchlist_total_count"] != len(CFF_2026_PRIORITY_WATCHLIST):
            print("FAIL cff_priority_radar self-check", file=sys.stderr)
            return 1
        print(
            f"PASS cff_priority_radar: watchlist={report['watchlist_total_count']} "
            f"(interventional_snapshot={report['interventional_snapshot_count']} "
            f"within raw={report['snapshot_raw_interventional_count']}/prisma={report['prisma_included_interventional_count']}; "
            f"observational_or_external={report['observational_or_external_count']}); "
            f"falsify_controls={len(f_rep['negative_controls'])}/{len(f_rep['negative_controls'])}"
        )
        return 0

    if args.json:
        print(json.dumps(report, indent=2, ensure_ascii=False))
    else:
        print(f"=== {report['tool']} ({report['epistemic_tier']}) ===")
        print(f"Boundary: {report['clinical_boundary']}")
        print(f"Disclaimer: {report['nlm_disclaimer']}")
        print(
            f"Interventional snapshot watchlist: {report['interventional_snapshot_count']} "
            f"| External/Observational: {report['observational_or_external_count']}"
        )
        for r in report["interventional_snapshot_watchlist"]:
            print(f"  - [{r['nct_id']}] ({r['priority_pillar']}, {r['overall_status']}, prisma_940={r['included_in_prisma_940']}) {r['program_label']}")
        for r in report["observational_or_external_watchlist"]:
            print(f"  - [{r['nct_id']}] ({r['priority_pillar']}, {r['study_type']}, {r['overall_status']}) {r['program_label']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
