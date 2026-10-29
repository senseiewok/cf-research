#!/usr/bin/env python3
"""Comprehensive Unit Tests for Airway Basal Stem Cell ITGA6 Tropism & Clonal Regeneration Model (T-0498).

Verifies:
1. Receptor binding Michaelis-Menten kinetics and tight-junction permeability opening.
2. Targeting index selectivity (ITGA6-targeted vs non-targeted delivery).
3. Three-compartment clonal renewal ODE integration across 365+ days.
4. Clonal stem persistence vs transient non-stem shedding decay.
5. Primary evidence attribution, NLM disclaimer, and research boundaries.
6. Negative controls (zero vector dose, closed tight-junction window).
7. Epistemic robustness under +/-25% parameter perturbation (AH-11).
"""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from basal_stem_tropism_persistence import (
    CURATIVE_THRESHOLD_PCT_WT,
    DEFAULT_KD_ITGA6_TARGETED_NM,
    DEFAULT_KD_NONTARGETED_NM,
    DEFAULT_PO_WT,
    DOI_COONEY_BASAL_LENTIVIRAL,
    DOI_IVLEVA_BASAL_EV,
    DOI_LIANG_PRIME_EDITING,
    DOI_TATA_EPITHELIAL_PLASTICITY,
    DOI_VIJAYKUMAR_BASAL_TRANSPLANT,
    NLM_ATTRIBUTION,
    PDB_ACTIVE_CFTR,
    RESEARCH_BOUNDARY,
    calculate_basal_transduction_flux,
    calculate_tight_junction_permeability,
    evaluate_parameter_sensitivity,
    simulate_basal_tropism_persistence,
    solve_three_compartment_clonal_ode,
)


class TestBasalStemTropismPersistence(unittest.TestCase):
    """Test suite for T-0498 airway basal stem cell targeting and clonal regeneration model."""

    def test_attribution_and_primary_grounding(self) -> None:
        """Verify NLM attribution, research boundary, and primary evidence anchors."""
        self.assertIn("ClinicalTrials.gov does not evaluate safety or efficacy", NLM_ATTRIBUTION)
        self.assertIn("Preclinical computational biophysical simulation", RESEARCH_BOUNDARY)
        self.assertEqual(PDB_ACTIVE_CFTR, "5UAK")
        self.assertEqual(DOI_COONEY_BASAL_LENTIVIRAL, "10.1016/j.ymthe.2018.04.014")
        self.assertEqual(DOI_VIJAYKUMAR_BASAL_TRANSPLANT, "10.1172/jci.insight.206375")
        self.assertEqual(DOI_IVLEVA_BASAL_EV, "10.3389/fbioe.2026.1855907")
        self.assertEqual(DOI_TATA_EPITHELIAL_PLASTICITY, "10.1038/nature12777")
        self.assertEqual(DOI_LIANG_PRIME_EDITING, "10.1016/j.xcrm.2024.101544")

    def test_tight_junction_permeability_window(self) -> None:
        """Verify tight-junction opening function theta_TJ(w) satisfies physical bounds [0, 1]."""
        self.assertEqual(calculate_tight_junction_permeability(0.0), 0.0)
        self.assertAlmostEqual(calculate_tight_junction_permeability(1.0, k_tj=1.0), 0.5, places=4)
        self.assertGreater(calculate_tight_junction_permeability(4.0), calculate_tight_junction_permeability(1.0))
        self.assertLess(calculate_tight_junction_permeability(100.0), 1.0)

    def test_receptor_mediated_targeting_selectivity(self) -> None:
        """Verify ITGA6-targeted vectors exhibit superior basolateral transduction flux vs non-targeted."""
        targeted = calculate_basal_transduction_flux(vector_apical_nm=5.0, window_hours=2.0, is_itga6_targeted=True)
        nontargeted = calculate_basal_transduction_flux(vector_apical_nm=5.0, window_hours=2.0, is_itga6_targeted=False)

        self.assertTrue(targeted.is_targeted)
        self.assertFalse(nontargeted.is_targeted)
        self.assertLess(targeted.kd_nm, nontargeted.kd_nm)
        self.assertGreater(targeted.j_basal_transduction_rate, nontargeted.j_basal_transduction_rate)
        self.assertGreater(targeted.targeting_index_ratio, 10.0)

    def test_basal_transduction_dose_response(self) -> None:
        """Verify monotonic increase with dose and asymptotic saturation."""
        low_dose = calculate_basal_transduction_flux(vector_apical_nm=1.0, window_hours=2.0)
        mid_dose = calculate_basal_transduction_flux(vector_apical_nm=10.0, window_hours=2.0)
        high_dose = calculate_basal_transduction_flux(vector_apical_nm=100.0, window_hours=2.0)

        self.assertLess(low_dose.j_basal_transduction_rate, mid_dose.j_basal_transduction_rate)
        self.assertLess(mid_dose.j_basal_transduction_rate, high_dose.j_basal_transduction_rate)
        # Check asymptotic limit <= theta_TJ * J0
        self.assertLessEqual(high_dose.j_basal_transduction_rate, high_dose.theta_tj)

    def test_three_compartment_ode_integration(self) -> None:
        """Verify numerical stability and positive trajectory across 365 days."""
        trajectory = solve_three_compartment_clonal_ode(
            b_star_0=2.5e4,  # 10% basal transduction
            t_star_0=1.0e4,
            d_star_0=5.0e3,
            total_days=365.0,
            dt=0.25,
        )
        self.assertGreater(len(trajectory), 300)
        first_point = trajectory[0]
        last_point = trajectory[-1]

        self.assertEqual(first_point.day, 0.0)
        self.assertAlmostEqual(last_point.day, 365.0, delta=1.0)
        self.assertGreater(last_point.differentiated_transduced, 0.0)
        self.assertGreater(last_point.restored_isc_pct_wt, 0.0)

    def test_clonal_persistence_vs_non_stem_decay(self) -> None:
        """Verify basal stem targeting provides sustained expression whereas non-stem delivery decays."""
        stem_sim = simulate_basal_tropism_persistence(vector_dose_apical_nm=10.0, window_hours=2.0, is_itga6_targeted=True)
        nonstem_sim = simulate_basal_tropism_persistence(vector_dose_apical_nm=10.0, window_hours=2.0, is_itga6_targeted=False)

        # Stem targeting maintains robust clonal renewal expansion (Day 365 / Day 30 > 1.0)
        self.assertGreater(stem_sim.clonal_persistence_index_365_over_30, 1.0)
        # Non-stem delivery suffers significant loss due to epithelial shedding turnover
        self.assertLess(nonstem_sim.clonal_persistence_index_365_over_30, 0.40)
        # Absolute current at day 365 is much higher with stem targeting
        self.assertGreater(stem_sim.day_365_isc_pct_wt, nonstem_sim.day_365_isc_pct_wt * 5.0)

    def test_curative_threshold_achievement(self) -> None:
        """Verify targeted basal stem delivery exceeds >= 10% WT CFTR curative threshold at day 365."""
        res = simulate_basal_tropism_persistence(vector_dose_apical_nm=15.0, window_hours=2.5, is_itga6_targeted=True)
        self.assertTrue(res.symptom_relief_threshold_met)
        self.assertGreaterEqual(res.day_365_isc_pct_wt, CURATIVE_THRESHOLD_PCT_WT)
        self.assertGreaterEqual(res.curative_duration_days, 250.0)

    def test_negative_control_zero_dose(self) -> None:
        """Negative control: zero vector dose produces 0% CFTR restoration."""
        null_res = simulate_basal_tropism_persistence(vector_dose_apical_nm=0.0, window_hours=2.0)
        self.assertEqual(null_res.day_365_isc_pct_wt, 0.0)
        self.assertEqual(null_res.initial_basal_transduced_fraction, 0.0)
        self.assertTrue(null_res.negative_control_passed)
        self.assertFalse(null_res.symptom_relief_threshold_met)

    def test_negative_control_closed_tight_junctions(self) -> None:
        """Negative control: intact unpermeabilized tight junctions (window = 0) blocks basal transduction."""
        closed_res = simulate_basal_tropism_persistence(vector_dose_apical_nm=20.0, window_hours=0.0)
        self.assertEqual(closed_res.day_365_isc_pct_wt, 0.0)
        self.assertTrue(closed_res.negative_control_passed)

    def test_dedifferentiation_contribution(self) -> None:
        """Verify dedifferentiation per Tata et al. elevates long-term clonal basal replenishment."""
        with_dediff = simulate_basal_tropism_persistence(vector_dose_apical_nm=10.0, enable_dedifferentiation=True)
        without_dediff = simulate_basal_tropism_persistence(vector_dose_apical_nm=10.0, enable_dedifferentiation=False)

        self.assertTrue(with_dediff.dedifferentiation_supported)
        self.assertFalse(without_dediff.dedifferentiation_supported)
        self.assertGreaterEqual(with_dediff.day_365_isc_pct_wt, without_dediff.day_365_isc_pct_wt)

    def test_parametric_sensitivity_ah11(self) -> None:
        """Verify AH-11 epistemic robustness under +/-25% parameter perturbation."""
        sens = evaluate_parameter_sensitivity(base_dose=10.0, base_window=2.0, perturbation_pct=0.25)
        self.assertEqual(sens["rule"], "AH-11")
        self.assertEqual(sens["robustness_status"], "PASS")
        self.assertTrue(sens["dose_sensitivity"]["monotonic"])
        self.assertTrue(sens["tight_junction_window_sensitivity"]["monotonic"])

    def test_cli_execution_json(self) -> None:
        """Verify CLI execution with --json returns valid JSON and zero exit code."""
        cmd = [sys.executable, str(HERE / "basal_stem_tropism_persistence.py"), "--dose", "8.0", "--window", "1.5", "--json"]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, f"CLI stderr: {proc.stderr}")
        data = json.loads(proc.stdout)
        self.assertIn("day_365_isc_pct_wt", data)
        self.assertIn("clonal_persistence_index_365_over_30", data)
        self.assertTrue(data["negative_control_passed"])


if __name__ == "__main__":
    unittest.main()
