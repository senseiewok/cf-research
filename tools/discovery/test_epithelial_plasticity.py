#!/usr/bin/env python3
"""Unit tests for epithelial_plasticity.py (T-0178)."""
import unittest
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import epithelial_plasticity as ep


class TestEpithelialPlasticityModel(unittest.TestCase):
    def test_baseline_parameters(self):
        self.assertIn("basal_stem_cells", ep.AIRWAY_EPITHELIAL_BASELINE["cell_fractions"])
        self.assertIn("secretory_club_cells", ep.AIRWAY_EPITHELIAL_BASELINE["cell_fractions"])
        self.assertEqual(ep.AIRWAY_EPITHELIAL_BASELINE["therapeutic_cftr_threshold_pct"], 15.0)

    def test_dedifferentiation_success(self):
        res = ep.simulate_dedifferentiation_bypass(
            apical_editing_efficiency=0.40,
            notch_inhibitor_pulse_days=5.0,
            total_days=180.0
        )
        self.assertTrue(res["meets_therapeutic_threshold"])
        self.assertGreaterEqual(res["final_basal_stem_cell_chimerism_pct"], 15.0)
        self.assertEqual(res["verdict"], "PERMANENT_CURATIVE_REPOPULATION")
        # Long-term persistence at day 180 should stay well above zero
        self.assertGreater(res["trajectory_sample_day180"], 15.0)

    def test_zero_notch_pulse_leads_to_ephemeral_expression(self):
        # Without Notch inhibitor pulse, no dedifferentiation occurs: edited apical cells turn over and wash out
        res = ep.simulate_dedifferentiation_bypass(
            apical_editing_efficiency=0.40,
            notch_inhibitor_pulse_days=0.0,
            total_days=180.0
        )
        self.assertEqual(res["final_basal_stem_cell_chimerism_pct"], 0.0)
        self.assertFalse(res["meets_therapeutic_threshold"])
        self.assertEqual(res["verdict"], "EPHEMERAL_EXPRESSION_NO_REPOPULATION")

    def test_invalid_parameters(self):
        with self.assertRaises(ValueError):
            ep.simulate_dedifferentiation_bypass(apical_editing_efficiency=1.5)
        with self.assertRaises(ValueError):
            ep.simulate_dedifferentiation_bypass(notch_inhibitor_pulse_days=-1.0)


if __name__ == "__main__":
    unittest.main()
