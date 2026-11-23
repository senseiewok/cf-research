#!/usr/bin/env python3
"""Unit and falsification tests for tools/discovery/mucus_barrier_screener.py (T-0921)."""
from __future__ import annotations

import json
import re
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import mucus_barrier_screener as mbs  # noqa: E402


class TestMucusBarrierScreener(unittest.TestCase):
    def test_1_dense_peg_brush_mpp_crosses_cf_sputum_before_mcc(self) -> None:
        res = mbs.evaluate_mucus_barrier(
            mbs.VectorBiophysicalInput(
                diameter_nm=80.0,
                zeta_mv=-2.0,
                peg_density_chains_nm2=0.28,
                peg_mw_kda=2.0,
                pore_diameter_nm=140.0,
                mucus_thickness_um=25.0,
            )
        )
        cm = res["computed_metrics"]
        self.assertEqual(cm["peg_rf_over_d"]["regime"], "dense_brush")
        self.assertTrue(cm["crosses_before_mcc_window"]["value"])
        self.assertEqual(cm["crosses_before_mcc_window"]["verdict"], "FEASIBLE_MUCUS_PENETRATING_CANDIDATE")
        self.assertGreaterEqual(len(res["not_modeled_in_vivo"]), 4)

    def test_2_oversized_or_cationic_vectors_are_trapped(self) -> None:
        oversized = mbs.evaluate_mucus_barrier(
            mbs.VectorBiophysicalInput(diameter_nm=220.0, pore_diameter_nm=140.0)
        )
        self.assertFalse(oversized["computed_metrics"]["crosses_before_mcc_window"]["value"])
        self.assertEqual(
            oversized["computed_metrics"]["crosses_before_mcc_window"]["verdict"],
            "STERICALLY_TRAPPED_BY_SPUTUM_MESH",
        )

        cationic = mbs.evaluate_mucus_barrier(
            mbs.VectorBiophysicalInput(diameter_nm=90.0, zeta_mv=32.0, peg_density_chains_nm2=0.0)
        )
        self.assertFalse(cationic["computed_metrics"]["crosses_before_mcc_window"]["value"])
        self.assertEqual(
            cationic["computed_metrics"]["crosses_before_mcc_window"]["verdict"],
            "ELECTROSTATICALLY_TRAPPED_ON_MUCIN_FIBERS",
        )

    def test_3_falsification_suite_passes_all_negative_controls(self) -> None:
        rep = mbs.run_falsification_suite()
        self.assertTrue(rep["falsify_passed"], msg=str(rep))
        self.assertEqual(len(rep["negative_controls"]), 3)

    def test_4_empirical_mpt_benchmark_matches_published_regimes(self) -> None:
        b_rep = mbs.evaluate_empirical_mpt_benchmark()
        self.assertEqual(b_rep["total_cases"], 5)
        self.assertEqual(b_rep["passed_cases"], 5, msg=str(b_rep["cases"]))
        self.assertTrue(b_rep["passed"])

    def test_5_dossier_table_matches_screener_outputs(self) -> None:
        dossier_path = HERE.parents[1] / "landscape" / "airway-mucus-barrier-biophysics-2026.md"
        text = dossier_path.read_text(encoding="utf-8")
        self.assertNotIn("Stokes-Einstein-Ogston-Brenner", text)
        self.assertNotIn("SEOB", text)

        rows_by_id: dict[str, list[str]] = {}
        for line in text.splitlines():
            stripped = line.strip()
            if not stripped.startswith("| `"):
                continue
            cells = [c.strip().strip("`") for c in stripped.strip("|").split("|")]
            if len(cells) >= 10:
                rows_by_id[cells[0]] = cells

        self.assertEqual(len(rows_by_id), len(mbs.EMPIRICAL_MPT_BENCHMARK_CASES))
        for case in mbs.EMPIRICAL_MPT_BENCHMARK_CASES:
            case_id = case["case_id"]
            self.assertIn(case_id, rows_by_id)
            cells = rows_by_id[case_id]
            inp = case["input"]
            pkt = mbs.evaluate_mucus_barrier(inp)
            cm = pkt["computed_metrics"]

            self.assertAlmostEqual(float(cells[2]), inp.diameter_nm, places=6)
            self.assertAlmostEqual(float(cells[3]), inp.zeta_mv, places=6)
            self.assertAlmostEqual(float(cells[4]), inp.pore_diameter_nm, places=6)
            self.assertAlmostEqual(float(cells[5]), inp.mucus_thickness_um, places=6)
            self.assertAlmostEqual(
                float(cells[6]),
                float(cm["d_eff_over_d_water"]["value"]),
                places=6,
            )
            transit_expected = cm["mean_transit_time_s"]["value"]
            if transit_expected is None:
                self.assertIn("sterically blocked", cells[8].lower())
            else:
                self.assertAlmostEqual(float(cells[8]), float(transit_expected), places=2)
            self.assertEqual(cells[9], str(cm["crosses_before_mcc_window"]["verdict"]))

        # Verify companion evidence corpus and claims ledger match screener outputs too.
        landscape_dir = HERE.parents[1] / "landscape"
        evidence_text = (landscape_dir / "airway-mucus-barrier-biophysics-2026-evidence.txt").read_text(encoding="utf-8")
        claims_text = (landscape_dir / "airway-mucus-barrier-biophysics-2026-claims.json").read_text(encoding="utf-8")
        stale_tokens = ("Stokes-Einstein-Ogston-Brenner", "SEOB", "0.036923", "0.000018", "0.000516")
        for label, body in (("markdown", text), ("evidence", evidence_text), ("claims", claims_text)):
            for token in stale_tokens:
                self.assertNotIn(token, body, msg=f"Stale token {token!r} found in {label}")

        evidence_lines = {}
        for line in evidence_text.splitlines():
            m = re.match(r"^Regime\s+\d+\s+\(([^)]+)\):\s+(.+)$", line.strip())
            if m:
                evidence_lines[m.group(1)] = m.group(2)
        self.assertEqual(set(evidence_lines), {c["case_id"] for c in mbs.EMPIRICAL_MPT_BENCHMARK_CASES})

        claims_by_id = {entry["id"]: entry["claim"] for entry in json.loads(claims_text)}
        claim_ids = (
            "mpt-regime-1-80nm-peg",
            "mpt-regime-2-100nm-hybrid",
            "mpt-regime-3-25nm-aav",
            "mpt-regime-4-cationic-trapped",
            "mpt-regime-5-steric-ratio-computed",
        )
        for case, claim_id in zip(mbs.EMPIRICAL_MPT_BENCHMARK_CASES, claim_ids):
            case_id = case["case_id"]
            pkt = mbs.evaluate_mucus_barrier(case["input"])
            cm = pkt["computed_metrics"]
            expected_d = float(cm["d_eff_over_d_water"]["value"])
            expected_verdict = str(cm["crosses_before_mcc_window"]["verdict"])

            ev_m = re.search(r"D_eff/D_w of (\d+(?:\.\d+)?).*classified as ([A-Z_]+)\.$", evidence_lines[case_id])
            self.assertIsNotNone(ev_m, msg=f"Unparsed evidence line for {case_id}")
            self.assertAlmostEqual(float(ev_m.group(1)), expected_d, places=6)
            self.assertEqual(ev_m.group(2), expected_verdict)

            cl_m = re.search(r"D_eff/D_w of (\d+(?:\.\d+)?)", claims_by_id[claim_id])
            self.assertIsNotNone(cl_m, msg=f"Missing D_eff/D_w in claim {claim_id}")
            self.assertAlmostEqual(float(cl_m.group(1)), expected_d, places=6)


if __name__ == "__main__":
    unittest.main()
