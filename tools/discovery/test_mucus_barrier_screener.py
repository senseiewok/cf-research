#!/usr/bin/env python3
"""Unit and falsification tests for tools/discovery/mucus_barrier_screener.py (T-0921)."""
from __future__ import annotations

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


if __name__ == "__main__":
    unittest.main()
