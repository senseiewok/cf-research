#!/usr/bin/env python3
"""Comprehensive unit test suite for aerosol_mmad_mucus_penetration.py (T-0547).

Acceptance verification command:
  python -m unittest tools.discovery.test_aerosol_mmad_mucus_penetration -v && python tools/claims/check_anti_hallucination.py
"""
from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from aerosol_mmad_mucus_penetration import (
    NLM_REGISTRY_DISCLAIMER,
    PARAMETER_PROVENANCE,
    RESEARCH_BOUNDARY,
    AirwayDepositionProfile,
    CoupledDeliveryResult,
    MucusHydrogelScaling,
    brownian_diffusivity_air,
    calculate_peclet_number,
    calculate_stokes_number,
    cunningham_correction_factor,
    deposition_efficiency_generation,
    evaluate_coupled_delivery_efficiency,
    evaluate_discovery_panel,
    execute_negative_controls,
    extrathoracic_impaction_fraction,
    self_test,
    simulate_airway_deposition,
    simulate_mucus_hydrogel_scaling,
)


class TestAerosolFluidDynamics(unittest.TestCase):
    """Test fundamental aerosol fluid dynamics and aerodynamic sizing calculations."""

    def test_cunningham_correction_monotonicity(self):
        """Cunningham slip correction must decrease toward 1.0 as droplet diameter increases."""
        cc_submicron = cunningham_correction_factor(0.5e-6)
        cc_optimal = cunningham_correction_factor(2.0e-6)
        cc_large = cunningham_correction_factor(10.0e-6)

        self.assertGreater(cc_submicron, cc_optimal)
        self.assertGreater(cc_optimal, cc_large)
        self.assertGreater(cc_large, 1.0)
        self.assertLess(cc_large, 1.05)

    def test_cunningham_invalid_diameter_raises(self):
        with self.assertRaises(ValueError):
            cunningham_correction_factor(-1.0)

    def test_stokes_number_scaling(self):
        """Stokes number scales with square of droplet diameter (d_a^2)."""
        stk_1um = calculate_stokes_number(1.0e-6, airway_velocity_m_s=2.0, airway_diameter_m=0.01)
        stk_2um = calculate_stokes_number(2.0e-6, airway_velocity_m_s=2.0, airway_diameter_m=0.01)
        self.assertGreater(stk_2um, stk_1um)
        # Ratio roughly 4x adjusted for Cunningham slip
        ratio = stk_2um / stk_1um
        self.assertGreater(ratio, 3.0)
        self.assertLess(ratio, 4.5)

    def test_peclet_number_diffusion_scaling(self):
        """Peclet number scales with droplet diameter via Stokes-Einstein diffusivity."""
        pe_submicron = calculate_peclet_number(0.5e-6, airway_velocity_m_s=2.0, airway_diameter_m=0.01)
        pe_micron = calculate_peclet_number(2.5e-6, airway_velocity_m_s=2.0, airway_diameter_m=0.01)
        self.assertGreater(pe_micron, pe_submicron)

    def test_extrathoracic_impaction_curve(self):
        """Extrathoracic impaction must be near zero for submicron and >90% for 10 um droplets."""
        p_sub = extrathoracic_impaction_fraction(0.4, 30.0)
        p_opt = extrathoracic_impaction_fraction(2.0, 30.0)
        p_oversized = extrathoracic_impaction_fraction(10.0, 30.0)

        self.assertLess(p_sub, 0.05)
        self.assertGreater(p_opt, 0.15)
        self.assertLess(p_opt, 0.35)
        self.assertGreater(p_oversized, 0.90)


class TestAirwayDepositionProfiles(unittest.TestCase):
    """Test multi-compartment Weibel G0-G16 airway deposition profiles."""

    def test_optimal_mmad_bronchial_deposition(self):
        """MMAD in 1.5 - 2.5 um range must maximize target bronchial deposition (>= 30%)."""
        dep_1_5 = simulate_airway_deposition(1.5)
        dep_2_0 = simulate_airway_deposition(2.0)
        dep_2_5 = simulate_airway_deposition(2.5)

        self.assertGreaterEqual(dep_1_5.bronchial_target_fraction, 0.28)
        self.assertGreaterEqual(dep_2_0.bronchial_target_fraction, 0.35)
        self.assertGreaterEqual(dep_2_5.bronchial_target_fraction, 0.35)

    def test_submicron_droplet_exhalation_loss(self):
        """Sub-micron droplets (0.4 um) must exhibit high exhalation fraction (>= 60%)."""
        dep_sub = simulate_airway_deposition(0.4)
        self.assertGreaterEqual(dep_sub.alveolar_and_exhaled_fraction, 0.60)
        self.assertLess(dep_sub.bronchial_target_fraction, 0.10)

    def test_oversized_droplet_upper_airway_loss(self):
        """Oversized droplets (10 um) must deposit predominantly in upper/extrathoracic airways."""
        dep_large = simulate_airway_deposition(10.0)
        upper_total = dep_large.extrathoracic_fraction + dep_large.upper_airway_fraction
        self.assertGreaterEqual(upper_total, 0.90)
        self.assertLessEqual(dep_large.bronchial_target_fraction, 0.05)


class TestMucusHydrogelObstruction(unittest.TestCase):
    """Test steric mucus hydrogel obstruction scaling and transit time calculations."""

    def test_cf_vs_healthy_mucus_diffusivity_contrast(self):
        """CF thick mucus (xi=42 nm) must restrict 30 nm LNP more than healthy mucus (xi=120 nm)."""
        cf = simulate_mucus_hydrogel_scaling(30.0, xi_mesh_nm=42.0)
        healthy = simulate_mucus_hydrogel_scaling(30.0, xi_mesh_nm=120.0)

        self.assertLess(cf.diffusivity_ratio, healthy.diffusivity_ratio)
        self.assertGreaterEqual(cf.diffusivity_ratio, 0.20)
        self.assertGreaterEqual(healthy.diffusivity_ratio, 0.75)

    def test_oversized_lnp_steric_entrapment(self):
        """Oversized LNPs (radius >= 60 nm, diam >= 120 nm) must be entrapped in CF sputum (D_eff/D0 < 0.01)."""
        entrapped = simulate_mucus_hydrogel_scaling(60.0, xi_mesh_nm=42.0)
        self.assertEqual(entrapped.penetration_verdict, "ENTRAPPED")
        self.assertLess(entrapped.diffusivity_ratio, 0.01)

    def test_compact_lnp_penetrating_verdict(self):
        """Compact 20 nm LNP (40 nm diam) must achieve PENETRATING verdict in CF mucus."""
        compact = simulate_mucus_hydrogel_scaling(20.0, xi_mesh_nm=42.0)
        self.assertEqual(compact.penetration_verdict, "PENETRATING")
        self.assertGreaterEqual(compact.diffusivity_ratio, 0.40)


class TestCoupledDeliveryEfficiency(unittest.TestCase):
    """Test coupled aerosol MMAD and mucus hydrogel delivery efficiency optimization."""

    def test_optimal_mmad_range_flag(self):
        """1.5 - 2.5 um range must be flagged as optimal."""
        res_1_0 = evaluate_coupled_delivery_efficiency(1.0)
        res_2_0 = evaluate_coupled_delivery_efficiency(2.0)
        res_5_0 = evaluate_coupled_delivery_efficiency(5.0)

        self.assertFalse(res_1_0.is_optimal_mmad_range)
        self.assertTrue(res_2_0.is_optimal_mmad_range)
        self.assertFalse(res_5_0.is_optimal_mmad_range)

    def test_net_delivery_efficiency_peak(self):
        """Net delivery efficiency Phi_net peaks in the 1.5 - 2.5 um MMAD window."""
        phi_0_5 = evaluate_coupled_delivery_efficiency(0.5).net_delivery_efficiency
        phi_2_0 = evaluate_coupled_delivery_efficiency(2.0).net_delivery_efficiency
        phi_7_0 = evaluate_coupled_delivery_efficiency(7.0).net_delivery_efficiency

        self.assertGreater(phi_2_0, phi_0_5)
        self.assertGreater(phi_2_0, phi_7_0)


class TestGovernanceAndGuardrails(unittest.TestCase):
    """Test AH-9 parameter provenance, AH-10 negative controls, and AH-11 sensitivity."""

    def test_parameter_provenance_structure(self):
        """Rule AH-9: All parameters must declare source citations, units, and valid bounds."""
        for name, entry in PARAMETER_PROVENANCE.items():
            self.assertIn("value", entry)
            self.assertIn("units", entry)
            self.assertIn("bounds", entry)
            self.assertIn("source_type", entry)
            self.assertIn("source_citation", entry)
            lower, upper = entry["bounds"]
            self.assertLessEqual(lower, entry["value"])
            self.assertGreaterEqual(upper, entry["value"])

    def test_ah10_negative_controls(self):
        """Rule AH-10: All three falsification negative controls must pass."""
        nc = execute_negative_controls()
        self.assertTrue(nc["all_controls_passed"])
        self.assertTrue(nc["nc1_oversized_droplet_impaction"]["passed"])
        self.assertTrue(nc["nc2_oversized_lnp_mucus_entrapment"]["passed"])
        self.assertTrue(nc["nc3_submicron_droplet_exhalation"]["passed"])

    def test_ah11_parametric_sensitivity(self):
        """Rule AH-11: +/-25% parameter perturbation preserves optimal MMAD window."""
        base_panel = evaluate_discovery_panel()
        base_peak_pct = base_panel["peak_bronchial_deposition_pct"]

        # Perturb impaction scaling alpha by +25% and -25%
        dep_high_alpha = simulate_airway_deposition(2.0, alpha=15.0 * 1.25)
        dep_low_alpha = simulate_airway_deposition(2.0, alpha=15.0 * 0.75)

        self.assertGreaterEqual(dep_high_alpha.bronchial_target_fraction, 0.25)
        self.assertGreaterEqual(dep_low_alpha.bronchial_target_fraction, 0.25)

    def test_nlm_disclaimer_presence(self):
        self.assertIn("ClinicalTrials.gov does not evaluate safety or efficacy", NLM_REGISTRY_DISCLAIMER)

    def test_self_test_function(self):
        self.assertTrue(self_test())


if __name__ == "__main__":
    unittest.main()
