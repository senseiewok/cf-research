"""Unit tests for Endosomal Escape Kinetics & Proton Sponge Buffering Capacity (T-0557).

Primary Biophysical Grounding & Evidence:
  - DOI:10.1038/s41565-020-0669-6 (Selective organ targeting SORT lipid nanoparticles for tissue-specific mRNA delivery)
  - DOI:10.1038/s41587-023-01861-1 (Inhaled delivery of lipid nanoparticles to respiratory epithelia)
  - NCT05668741 (Phase 1/2 VX-522 inhaled mRNA LNP study in CF)
"""
from __future__ import annotations

import unittest
from tools.discovery.endosomal_escape_kinetics import (
    EndosomalEscapeParameters,
    compute_buffering_capacity,
    evaluate_endosomal_escape_kinetics,
    falsify_endosomal_escape_kinetics,
    simulate_endosomal_escape,
)


class EndosomalEscapeKineticsTests(unittest.TestCase):
    """Test suite verifying endosomal acidification, proton sponge buffering, and cytosolic escape."""

    def setUp(self) -> None:
        self.report = evaluate_endosomal_escape_kinetics()

    def test_report_structure_and_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0557")
        self.assertTrue(self.report["endosomal_escape_verified"])
        self.assertIn("SORT_Airway_LNP", self.report["formulations"])
        self.assertIn("Unbuffered_Cationic_LNP", self.report["formulations"])

    def test_buffering_capacity_computation(self) -> None:
        # Buffering capacity peaks near pH = pKa
        pka = 6.45
        c_lipid_mm = 35.0  # 35 mM
        h_at_pka = 10.0 ** (-pka)
        beta_peak = compute_buffering_capacity(h_at_pka, pka, c_lipid_mm)
        self.assertGreater(beta_peak, 15.0)

        # Buffering drops away from pKa (e.g. pH 4.0 or pH 8.0)
        h_low = 10.0 ** (-4.0)
        beta_low = compute_buffering_capacity(h_low, pka, c_lipid_mm)
        self.assertLess(beta_low, beta_peak)

        # Zero lipid concentration yields zero buffering
        self.assertEqual(compute_buffering_capacity(h_at_pka, pka, 0.0), 0.0)

    def test_formulation_escape_differential(self) -> None:
        sort_data = self.report["formulations"]["SORT_Airway_LNP"]
        cationic_data = self.report["formulations"]["Unbuffered_Cationic_LNP"]
        neutral_data = self.report["formulations"]["Neutral_Liposome"]

        # SORT Airway LNP achieves substantial escape yield (>= 8%)
        self.assertGreaterEqual(sort_data["cytosolic_escape_yield_pct"], 8.0)
        self.assertEqual(sort_data["escape_verdict"], "EFFICIENT_ENDOSOMAL_ESCAPE")
        self.assertGreater(sort_data["buffering_capacity_peak"], 15.0)

        # Unbuffered cationic LNP suffers predominantly lysosomal degradation
        self.assertLess(cationic_data["cytosolic_escape_yield_pct"], 3.0)
        self.assertEqual(cationic_data["escape_verdict"], "ENTRAPPED_LYSOSOMAL_DEGRADATION")

        # Neutral liposome has negligible escape
        self.assertLess(neutral_data["cytosolic_escape_yield_pct"], 1.0)

    def test_osmotic_pressure_and_membrane_rupture(self) -> None:
        sort_data = self.report["formulations"]["SORT_Airway_LNP"]
        # Counterion accumulation generates significant osmotic pressure (> 15 kPa)
        self.assertGreater(sort_data["peak_osmotic_pressure_kpa"], 15.0)
        self.assertGreater(sort_data["peak_rupture_probability"], 0.40)
        self.assertIsNotNone(sort_data["rupture_time_seconds"])

    def test_ah10_dynamic_falsification_negative_control(self) -> None:
        fals_result = falsify_endosomal_escape_kinetics()
        self.assertEqual(fals_result["rule"], "AH-10")
        self.assertEqual(fals_result["verdict"], "PASS")
        self.assertTrue(fals_result["negative_control_passed"])
        self.assertLess(fals_result["unbuffered_escape_yield_pct"], 1.0)


if __name__ == "__main__":
    unittest.main()
