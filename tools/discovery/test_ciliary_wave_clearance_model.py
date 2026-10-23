#!/usr/bin/env python3
"""Unit test suite for Airway Epithelial Ciliary Metachronal Wave Dynamics & Clearance Model (T-0315).

Acceptance Criteria:
  python -m unittest tools.discovery.test_ciliary_wave_clearance_model; python tools/discovery/ciliary_wave_clearance_model.py --check
"""
import math
import sys
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ciliary_wave_clearance_model import (
    BASELINE_HEALTHY_SOLIDS_WT_PCT,
    CRITICAL_MUCUS_SOLIDS_WT_PCT,
    CRITICAL_PCL_STASIS_HEIGHT_UM,
    HEALTHY_CBF_HZ,
    HEALTHY_PCL_HEIGHT_UM,
    NLM_REGISTRY_DISCLAIMER,
    calculate_effective_pcl_height,
    calculate_mucus_osmotic_pressure,
    calculate_mucus_viscosity,
    evaluate_cf_stasis,
    evaluate_stokeslet_flow,
    evaluate_two_phase_coupling,
    export_complete_ciliary_wave_manifest,
    run_checks,
)


class TestCiliaryWaveClearanceModel(unittest.TestCase):
    """Deterministic biophysical unit tests for T-0315."""

    def test_stokeslet_net_volume_transport(self):
        flow = evaluate_stokeslet_flow(HEALTHY_PCL_HEIGHT_UM, cbf_hz=HEALTHY_CBF_HZ)
        self.assertGreater(flow.forward_power_stroke_tip_velocity_um_s, 0.0)
        self.assertGreater(flow.volume_transport_rate_q_um2_s, 0.0)
        self.assertAlmostEqual(
            flow.mean_pcl_velocity_um_s * flow.pcl_height_um,
            flow.volume_transport_rate_q_um2_s,
            delta=1e-2,
        )

    def test_stokeslet_cbf_frequency_scaling(self):
        flow_8hz = evaluate_stokeslet_flow(cbf_hz=8.0)
        flow_12hz = evaluate_stokeslet_flow(cbf_hz=12.0)
        self.assertGreater(flow_12hz.forward_power_stroke_tip_velocity_um_s, flow_8hz.forward_power_stroke_tip_velocity_um_s)
        self.assertAlmostEqual(
            flow_12hz.volume_transport_rate_q_um2_s / flow_8hz.volume_transport_rate_q_um2_s,
            12.0 / 8.0,
            delta=0.01,
        )

    def test_mucus_viscosity_solids_power_law(self):
        eta_2_5 = calculate_mucus_viscosity(2.5)
        eta_5_0 = calculate_mucus_viscosity(5.0)
        eta_10_0 = calculate_mucus_viscosity(10.0)
        self.assertLess(eta_2_5, eta_5_0)
        self.assertLess(eta_5_0, eta_10_0)
        # Power law scaling should produce > 5x viscosity increase between 2.5% and 5.0%
        self.assertGreater(eta_5_0 / eta_2_5, 5.0)

    def test_mucus_osmotic_pressure_brush_gel_balance(self):
        pi_2_5 = calculate_mucus_osmotic_pressure(2.5)
        pi_8_0 = calculate_mucus_osmotic_pressure(8.0)
        # At healthy solids (2.5%), Pi_mucus < PCL brush resistance (~280 Pa)
        self.assertLess(pi_2_5, 280.0)
        # At critical solids (8.0%), Pi_mucus exceeds PCL resistance
        self.assertGreater(pi_8_0, 280.0)

    def test_pcl_compression_under_hyperconcentrated_mucus(self):
        h_healthy, _, _ = calculate_effective_pcl_height(2.5)
        h_cf, _, _ = calculate_effective_pcl_height(9.0)
        self.assertEqual(h_healthy, HEALTHY_PCL_HEIGHT_UM)
        self.assertLess(h_cf, CRITICAL_PCL_STASIS_HEIGHT_UM)

    def test_two_phase_coupling_healthy_clearance_velocity(self):
        res = evaluate_two_phase_coupling(BASELINE_HEALTHY_SOLIDS_WT_PCT)
        # Healthy physiological velocity in mm/min is 3 - 7 mm/min
        self.assertGreaterEqual(res.mucus_transport_velocity_mm_min, 3.0)
        self.assertLessEqual(res.mucus_transport_velocity_mm_min, 7.0)
        self.assertEqual(res.clearance_regime, "PHYSIOLOGICAL_RAPID_CLEARANCE")

    def test_cf_mucociliary_stasis_arrest(self):
        stasis_res = evaluate_cf_stasis(CRITICAL_MUCUS_SOLIDS_WT_PCT + 0.5)
        self.assertTrue(stasis_res.is_stasis)
        self.assertEqual(stasis_res.transport_velocity_mm_min, 0.0)
        self.assertEqual(stasis_res.effective_cbf_hz, 0.0)
        self.assertEqual(stasis_res.ciliary_coordination_factor, 0.0)

    def test_cf_subthreshold_partial_clearance(self):
        sub_res = evaluate_cf_stasis(4.5)
        self.assertFalse(sub_res.is_stasis)
        self.assertGreater(sub_res.transport_velocity_mm_min, 0.5)
        self.assertLess(sub_res.transport_velocity_mm_min, 6.0)

    def test_manifest_grounding_and_dois(self):
        manifest = export_complete_ciliary_wave_manifest()
        self.assertEqual(manifest["model_id"], "T-0315")
        self.assertEqual(manifest["answer_tier"], "T2")
        self.assertEqual(manifest["track"], "Track A")
        self.assertIn("NCT05668741", manifest["clinical_trial_grounding"])
        self.assertIn("10.1126/science.1223012", manifest["primary_literature_dois"])
        self.assertIn("10.1016/s0022-5193(86)80098-4", manifest["primary_literature_dois"])
        self.assertIn("10.1056/NEJMra1813799", manifest["primary_literature_dois"])
        self.assertIn("ClinicalTrials.gov does not evaluate safety or efficacy", manifest["nlm_registry_disclaimer"])

    def test_self_check_runner(self):
        code = run_checks()
        self.assertEqual(code, 0)


if __name__ == "__main__":
    unittest.main()
