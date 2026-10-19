#!/usr/bin/env python3
"""Unit tests for mechano_acoustic_delivery.py."""
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import mechano_acoustic_delivery as mad


class TestMechanoAcousticDelivery(unittest.TestCase):
    def test_yield_stress_monotonic_with_solids(self):
        tau_healthy = mad.compute_yield_stress_pa(2.0)
        tau_cf = mad.compute_yield_stress_pa(7.5)
        self.assertAlmostEqual(tau_healthy, 0.12, places=3)
        self.assertGreater(tau_cf, 2.5)
        with self.assertRaises(ValueError):
            mad.compute_yield_stress_pa(0.2)

    def test_thixotropic_conditioning_breakdown_and_window(self):
        static = mad.evaluate_thixotropic_conditioning(mad.BENCHMARK_REGIMES[0], solids_pct=7.5)
        hfcwo = mad.evaluate_thixotropic_conditioning(mad.BENCHMARK_REGIMES[1], solids_pct=7.5)
        bubble = mad.evaluate_thixotropic_conditioning(mad.BENCHMARK_REGIMES[3], solids_pct=7.5)

        self.assertAlmostEqual(static.steady_structure_lambda, 1.0, places=3)
        self.assertEqual(static.thixotropic_window_minutes, 0.0)
        self.assertFalse(static.asl_hydration_restored)

        # HFCWO and Acoustic + Microbubble reduce lambda, lower G', widen pores, and open a dosing window
        self.assertLess(hfcwo.steady_structure_lambda, static.steady_structure_lambda)
        self.assertLess(bubble.steady_structure_lambda, hfcwo.steady_structure_lambda)
        self.assertGreater(hfcwo.conditioned_pore_diameter_nm, static.conditioned_pore_diameter_nm)
        self.assertGreater(bubble.conditioned_pore_diameter_nm, hfcwo.conditioned_pore_diameter_nm)
        self.assertGreater(hfcwo.thixotropic_window_minutes, 5.0)
        self.assertTrue(bubble.asl_hydration_restored)
        self.assertGreaterEqual(bubble.predicted_asl_height_um, mad.PHYSIOLOGICAL_ASL_TARGET_UM)

    def test_vector_penetration_synergy_gains(self):
        lnp = mad.FRONTIER_VECTORS[0]  # 72 nm SORT-LNP
        ev = mad.FRONTIER_VECTORS[1]   # 85 nm Extracellular Vesicle
        minidna = mad.FRONTIER_VECTORS[2]  # 42 nm Minimal DNA vector

        syn_lnp = mad.evaluate_vector_penetration_synergy(lnp, mad.BENCHMARK_REGIMES[3], solids_pct=7.5)
        syn_ev = mad.evaluate_vector_penetration_synergy(ev, mad.BENCHMARK_REGIMES[3], solids_pct=7.5)
        syn_dna = mad.evaluate_vector_penetration_synergy(minidna, mad.BENCHMARK_REGIMES[3], solids_pct=7.5)

        self.assertGreater(syn_lnp.diffusivity_fold_gain, 2.0)
        self.assertGreater(syn_lnp.conditioned_epithelial_arrival_pct, syn_lnp.static_epithelial_arrival_pct)
        self.assertGreater(syn_ev.conditioned_epithelial_arrival_pct, syn_ev.static_epithelial_arrival_pct)
        # Smaller 42 nm minimal DNA vector achieves higher arrival than 85 nm EV
        self.assertGreater(syn_dna.conditioned_epithelial_arrival_pct, syn_ev.conditioned_epithelial_arrival_pct)

    def test_full_synergy_matrix_determinism(self):
        r1 = mad.run_full_synergy_matrix(7.5)
        r2 = mad.run_full_synergy_matrix(7.5)
        self.assertEqual(r1, r2)
        self.assertEqual(len(r1["conditioning_states"]), 4)
        self.assertEqual(len(r1["vector_synergies"]), 20)


if __name__ == "__main__":
    unittest.main()
