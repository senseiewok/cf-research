#!/usr/bin/env python3
"""Post-registration ClinicalTrials.gov API v2 `contactsLocationsModule` acquisition & PII stripper (`PRED-0018`).

Retrieves `protocolSection.identificationModule` and `protocolSection.contactsLocationsModule` once for
all 71 locked NCT IDs via `cf-evidence-loop` `Client` (`sources/catalog.yaml` `clinicaltrials-gov` `access: api`,
`<= 1.0 req/s`), verifies `71/71` `Status.FOUND` and `nctId` identity, strips `centralContacts`,
`overallOfficials`, `locations[].contacts`, `facility`, `city`, `state`, `zip`, and `geoPoint` in memory
before disk persistence, and validates every sanitized record with `parse_study_locations_record`.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path
import sys
from typing import Any
import urllib.parse

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import early_phase_geographic_footprint as bmod  # noqa: E402


def build_pii_stripped_fixture_from_client(
    client: Any,
    pred0015_doc: dict[str, Any],
    status_found_sentinel: Any,
    evaluated_on_local: str | None = None,
) -> dict[str, Any]:
    """Fetch all 71 locked NCT records via `client.get`, strip PII/facility fields in memory, and validate."""
    pe15 = pred0015_doc["primary_estimand_early_phase"]
    sys_ncts = list(pe15["systemic_ncts"])
    inh_ncts = list(pe15["inhaled_ncts"])
    dual_set = set(pe15["dual_route_inhaled_ncts"])
    all_ncts = sorted(set(sys_ncts) | set(inh_ncts))
    if len(all_ncts) != 71:
        raise bmod.EarlyPhaseGeographicBenchmarkError(f"Expected 71 locked NCTs, got {len(all_ncts)}")

    sys_set = set(sys_ncts)
    studies_out: dict[str, dict[str, Any]] = {}

    for nct_id in all_ncts:
        resp = client.get("clinicaltrials-gov", f"studies/{nct_id}")
        parsed_url = urllib.parse.urlparse(str(getattr(resp, "url", "")))
        if parsed_url.scheme != "https" or parsed_url.hostname != "clinicaltrials.gov":
            raise bmod.EarlyPhaseGeographicBenchmarkError(f"Unexpected upstream URL for {nct_id}: {getattr(resp, 'url', None)!r}")
        if getattr(resp, "status", None) != status_found_sentinel or not isinstance(getattr(resp, "data", None), dict):
            raise bmod.EarlyPhaseGeographicBenchmarkError(
                f"{nct_id}: non-FOUND response ({getattr(resp, 'status', None)!r}) from ClinicalTrials.gov API v2"
            )
        payload = resp.data
        proto = payload.get("protocolSection")
        if not isinstance(proto, dict):
            raise bmod.EarlyPhaseGeographicBenchmarkError(f"{nct_id}: missing protocolSection")
        ident = proto.get("identificationModule")
        if not isinstance(ident, dict) or str(ident.get("nctId") or "") != nct_id:
            raise bmod.EarlyPhaseGeographicBenchmarkError(f"{nct_id}: identificationModule.nctId mismatch")

        raw_cl_mod = proto.get("contactsLocationsModule")
        if raw_cl_mod is not None and not isinstance(raw_cl_mod, dict):
            raise bmod.EarlyPhaseGeographicBenchmarkError(f"{nct_id}: contactsLocationsModule must be dict or None")
        cl_mod_present = isinstance(raw_cl_mod, dict)
        cl_mod = raw_cl_mod if isinstance(raw_cl_mod, dict) else {}
        locs_key_present = cl_mod_present and ("locations" in cl_mod) and (cl_mod.get("locations") is not None)
        raw_locs = cl_mod.get("locations") if locs_key_present else None

        sanitized_locs: list[dict[str, str | None]] | None
        if raw_locs is None:
            sanitized_locs = None
        else:
            if not isinstance(raw_locs, list):
                raise bmod.EarlyPhaseGeographicBenchmarkError(f"{nct_id}: locations is not a list")
            sanitized_locs = []
            for loc_entry in raw_locs:
                if not isinstance(loc_entry, dict):
                    raise bmod.EarlyPhaseGeographicBenchmarkError(f"{nct_id}: location entry is not a dict")
                raw_country = loc_entry.get("country")
                raw_status = loc_entry.get("status")
                clean_country = bmod.sanitize_registry_text(raw_country) if isinstance(raw_country, str) else raw_country
                clean_status = bmod.sanitize_registry_text(raw_status) if isinstance(raw_status, str) else raw_status
                sanitized_locs.append({"country": clean_country, "status": clean_status})

        expected_arm = "systemic_cftr" if nct_id in sys_set else "inhaled_cftr"
        expected_dual = nct_id in dual_set
        rec_out: dict[str, Any] = {
            "nct_id": nct_id,
            "arm": expected_arm,
            "dual_route_modulator_co_admin": expected_dual,
            "contacts_locations_module_present": cl_mod_present,
            "locations_key_present": locs_key_present,
            "locations": sanitized_locs,
        }
        bmod.parse_study_locations_record(
            rec_out,
            nct_id=nct_id,
            expected_arm=expected_arm,
            expected_dual_route=expected_dual,
        )
        studies_out[nct_id] = rec_out

    if len(studies_out) != 71:
        raise bmod.EarlyPhaseGeographicBenchmarkError(f"Incomplete fetch: {len(studies_out)}/71")

    fixture_doc = {
        "schema_version": "pred0018-early-phase-locations-fixture-v1",
        "catalog_source_id": "clinicaltrials-gov",
        "api_endpoint_pattern": "https://clinicaltrials.gov/api/v2/studies/{nct_id}",
        "fetched_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "evaluated_on_local": evaluated_on_local or dt.date.today().isoformat(),
        "nlm_attribution": bmod.NLM_ATTRIBUTION,
        "research_boundary": bmod.RESEARCH_BOUNDARY,
        "quarantine_notice": (
            "External registry country and recruitment-status strings sanitized for bidi/default-ignorable controls "
            "and closing-tag lookalikes; centralContacts, overallOfficials, locations[].contacts, facility, city, "
            "state, zip, and geoPoint stripped in memory prior to disk write."
        ),
        "study_count": len(studies_out),
        "studies": studies_out,
    }
    bmod.assert_zero_at_sign_anywhere(fixture_doc, path_label="fixture_doc")
    return fixture_doc


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pred0015-results", type=Path, default=bmod.DEFAULT_PRED0015_RESULTS_JSON)
    ap.add_argument("--out", type=Path, default=bmod.DEFAULT_FIXTURE_JSON)
    ap.add_argument("--evidence-scripts", type=Path, default=None)
    args = ap.parse_args(argv)

    ev_dir = args.evidence_scripts or (
        REPO_ROOT.parent / "cf-skills" / ".claude" / "skills" / "cf-evidence-loop" / "scripts"
    )
    if str(ev_dir) not in sys.path:
        sys.path.insert(0, str(ev_dir))
    from evidence.http import Client  # type: ignore
    from evidence.record import Status  # type: ignore

    p15_doc = json.loads(args.pred0015_results.read_text(encoding="utf-8"))
    client = Client()
    doc = build_pii_stripped_fixture_from_client(client, p15_doc, Status.FOUND)
    args.out.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"WROTE {args.out} ({doc['study_count']} studies)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
