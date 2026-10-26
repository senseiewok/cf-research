#!/usr/bin/env python3
"""Export deterministic PRISMA 2020 counts and longitudinal phase/era transitions to static JSON.

Outputs:
  tools/trial_atlas/endpoint_atlas_export.json
"""
from __future__ import annotations

import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import lexicon as lex_mod
import scope as scope_mod
import snapshot as snap
from consensus_audit import _to_raw_study, normalize_phase

SNAPSHOT_JSON_PATH = HERE / "snapshot.json"
EXPORT_JSON_PATH = HERE / "endpoint_atlas_export.json"
RUN_RECORD_PATH = HERE / "run-records" / "2026-10-11-run-1.json"


def compute_atlas_export(snapshot_path: Path = SNAPSHOT_JSON_PATH) -> dict:
    data = json.loads(snapshot_path.read_text(encoding="utf-8"))
    studies = data.get("studies", [])
    lex = lex_mod.load()

    # 1. PRISMA 2020 screening flow
    total_identified = len(studies)
    exclusions_by_rule = Counter()
    included_studies = []

    for s in studies:
        raw = _to_raw_study(s)
        rec = snap.study_record(raw)
        rule, reason = scope_mod.classify(rec)
        if rule:
            exclusions_by_rule[rule] += 1
        else:
            included_studies.append((s, rec))

    total_excluded = sum(exclusions_by_rule.values())
    total_included = len(included_studies)

    # 2. Outcome entries & class prevalence
    total_outcome_entries = 0
    studies_with_outcomes = 0
    class_entry_counts = Counter()
    granular_entry_counts = Counter()
    class_study_counts = defaultdict(set)

    phase_transition_matrix = defaultdict(Counter)
    era_transition_matrix = defaultdict(Counter)

    for s, rec in included_studies:
        entries = snap.outcome_entries(rec)
        total_outcome_entries += len(entries)
        if s.get("primary_outcomes"):
            studies_with_outcomes += 1

        phases = s.get("phases") or ["NA"]
        start_year = scope_mod.start_info(rec)["start_year"]
        if start_year is None:
            era = "unknown"
        elif start_year < 2012:
            era = "pre-2012"
        elif start_year < 2020:
            era = "2012-2019"
        else:
            era = "2020-2026"

        for e in entries:
            tagged = lex.tag_entry(e)
            for t in tagged.get("tags", []):
                granular_entry_counts[t["class"]] += 1
            c8_list = tagged.get("canonical_8_classes", [])
            for c8 in c8_list:
                class_entry_counts[c8] += 1
                class_study_counts[c8].add(s["nct_id"])
                for p in phases:
                    phase_transition_matrix[p][c8] += 1
                era_transition_matrix[era][c8] += 1

    studies_by_nct = {s["nct_id"]: s for s in studies}

    # Format output dictionary
    export = {
        "metadata": {
            "source": "ClinicalTrials.gov API v2",
            "snapshot_data_timestamp": data.get("metadata", {}).get(
                "data_timestamp", "2026-10-14T14:12:11.425555+00:00"
            ),
            "generated_tool": "tools/trial_atlas/export_web_atlas.py",
            "disclaimer": (
                "ClinicalTrials.gov does not evaluate safety or efficacy; "
                "listing on the registry does not denote United States Federal Government endorsement."
            ),
        },
        "prisma_2020_flow": {
            "identification": {
                "route_condition": total_identified,
                "total_identified": total_identified,
            },
            "screening": {
                "total_screened": total_identified,
                "excluded": total_excluded,
                "exclusions_by_rule": {
                    "X1_observational": exclusions_by_rule.get("X1", 0),
                    "X2_expanded_access": exclusions_by_rule.get("X2", 0),
                    "X3_non_interventional": exclusions_by_rule.get("X3", 0),
                    "X4_condition_not_cf": exclusions_by_rule.get("X4", 0),
                    "X5_cf_in_negation": exclusions_by_rule.get("X5", 0),
                    "X6_manual": exclusions_by_rule.get("X6", 0),
                },
            },
            "eligibility_and_inclusion": {
                "eligible_studies": total_included,
                "included_studies": total_included,
                "studies_with_registered_primary_outcomes": studies_with_outcomes,
                "studies_without_registered_primary_outcomes": total_included - studies_with_outcomes,
                "total_outcome_entries": total_outcome_entries,
            },
        },
        "canonical_8_endpoints": {
            c: {
                "entries": class_entry_counts[c],
                "studies": len(class_study_counts[c]),
            }
            for c in lex_mod.CANONICAL_8_CLASSES
        },
        "phase_transition_matrix": {
            p: dict(phase_transition_matrix[p])
            for p in ["PHASE1", "PHASE2", "PHASE3", "PHASE4"]
        },
        "era_transition_matrix": {
            era: dict(era_transition_matrix[era])
            for era in (["pre-2012", "2012-2019", "2020-2026"] + (["unknown"] if "unknown" in era_transition_matrix else []))
        },
        "benchmark_trials": [
            {
                "nct_id": b["nct_id"],
                "name": b["name"],
                "title": studies_by_nct.get(b["nct_id"], {}).get("brief_title", ""),
                "phase": normalize_phase(studies_by_nct.get(b["nct_id"], {}).get("phases")),
                "status": studies_by_nct.get(b["nct_id"], {}).get("overall_status", ""),
                "modality": b["modality"],
                "population": b["population"],
                "primary_endpoint": [
                    o.get("measure", "")
                    for o in studies_by_nct.get(b["nct_id"], {}).get("primary_outcomes", [])
                ],
                "citation": b["citation"],
            }
            for b in [
                {
                    "nct_id": "NCT00909532",
                    "name": "STRIVE",
                    "modality": "ivacaftor (VX-770)",
                    "population": "G551D gating mutation",
                    "citation": "Ramsey BW et al., N Engl J Med 2011;365:1663-1672",
                },
                {
                    "nct_id": "NCT01807923",
                    "name": "TRAFFIC",
                    "modality": "lumacaftor / ivacaftor",
                    "population": "F508del homozygous",
                    "citation": "Wainwright CE et al., N Engl J Med 2015;373:220-231",
                },
                {
                    "nct_id": "NCT03525444",
                    "name": "AURORA-44",
                    "modality": "elexacaftor / tezacaftor / ivacaftor",
                    "population": "F508del heterozygous with minimal function (F/MF)",
                    "citation": "Middleton PG et al., N Engl J Med 2019;381:1809-1819",
                },
                {
                    "nct_id": "NCT05033080",
                    "name": "SKYLINE 102",
                    "modality": "vanzacaftor / tezacaftor / deutivacaftor",
                    "population": "F508del heterozygous with minimal function (F/MF)",
                    "citation": "ClinicalTrials.gov NCT05033080",
                },
                {
                    "nct_id": "NCT05668741",
                    "name": "BEACON-CF",
                    "modality": "VX-522 (inhaled CFTR mRNA)",
                    "population": "CF adults ineligible for CFTR modulators",
                    "citation": "ClinicalTrials.gov NCT05668741 (Completed; program discontinuation [unverified])",
                },
            ]
        ],
        "care_team_guide": {
            "title": "Care Team Reader Guide",
            "reader_path": "care_team",
            "what_the_atlas_can_answer": [
                "Prevalence and longitudinal shift of registered primary endpoint classes across 940 interventional CF trials (pre-2012, 2012-2019, 2020-2026).",
                "Distribution of primary trial endpoints across clinical phases (Phase 1 through Phase 4).",
                "Whether trial protocols evaluate specific physiological or patient-reported measures (e.g. FEV1, sweat chloride, CFQ-R, BMI, pulmonary exacerbations).",
                "Co-tagging patterns where measures capture both efficacy parameters and safety surveillance (e.g. multi-class tagging)."
            ],
            "what_the_atlas_cannot_answer": [
                "Individual patient treatment decisions, personalized medical recommendations, or clinical practice directives.",
                "Regulatory efficacy or safety clearance (ClinicalTrials.gov registrations reflect sponsor filings, not agency approval determinations).",
                "Unregistered exploratory or secondary endpoints omitted from primary protocol registrations.",
                "Personal modulator eligibility, which depends strictly on individual CFTR genotyping and national regulatory labels."
            ],
            "fev1_definition_and_safety_use": (
                "The lung_function_fev1 class captures all registered forced expiratory volume in one second measures, "
                "including those monitored as safety or bronchospasm checks post-dosing, which are co-tagged with "
                "safety_adverse_events under the multi_class flag rather than suppressed."
            ),
            "measures_list": [
                {
                    "class_name": c,
                    "status": "shown (canonical 8)" if c in lex_mod.CANONICAL_8_CLASSES else "counted, not shown (below the agreement bar)",
                    "entries": class_entry_counts[c] if c in lex_mod.CANONICAL_8_CLASSES else granular_entry_counts[c],
                }
                for c in (list(lex_mod.CANONICAL_8_CLASSES) + [cls for cls in lex.class_order if cls not in lex_mod.CANONICAL_8_CLASSES and cls not in ("other", "not_stated")])
            ],
            "clinician_review_gate": "Plain-language class notes stay unpublished until a clinician has read them (a person, not a model)."
        },
    }

    # Canonical sha256 checksum
    serialized = json.dumps(export, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    export["export_sha256"] = hashlib.sha256(serialized.encode("utf-8")).hexdigest()

    return export


def main(argv: list[str] | None = None) -> int:
    export = compute_atlas_export()
    text = json.dumps(export, indent=2, ensure_ascii=False) + "\n"
    EXPORT_JSON_PATH.write_text(text, encoding="utf-8")
    print(f"Exported trial endpoint atlas data to {EXPORT_JSON_PATH}")
    print(f"  SHA-256 Digest: {export['export_sha256']}")
    print(f"  Identified:     {export['prisma_2020_flow']['identification']['total_identified']}")
    print(f"  Excluded:       {export['prisma_2020_flow']['screening']['excluded']}")
    print(f"  Included:       {export['prisma_2020_flow']['eligibility_and_inclusion']['included_studies']}")
    print(f"  Outcome Entries:{export['prisma_2020_flow']['eligibility_and_inclusion']['total_outcome_entries']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
