#!/usr/bin/env python3
"""Unit tests for macrophage_efferocytosis_model.py (T-0251)."""

from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from tools.discovery.macrophage_efferocytosis_model import (
    DEFAULT_J_INFLUX_CELLS_ML_HR,
    DEFAULT_K_APO_PER_HR,
    DEFAULT_K_CAT_OVER_KM_NE_PER_NM_HR,
    DEFAULT_K_EFF_0_PER_MAC_HR,
    DEFAULT_K_NET_PER_HR,
    DEFAULT_K_SEC_NEC_PER_HR,
    DEFAULT_KM_EFF_CELLS_ML,
    DEFAULT_M_PHI_DENSITY_CELLS_ML,
    DEFAULT_MERTK_0_RECEPTORS_CELL,
    DEFAULT_NE_CONC_NM,
    DEFAULT_PH_CF,
    DEFAULT_PH_WT,
    DEFAULT_PK_LYS,
    NLM_REGISTRY_DISCLAIMER,
    RESEARCH_BOUNDARY,
    EfferocytosisParams,
    calculate_efferocytosis_rate,
    calculate_mertk_steady_state,
    calculate_phago_efficiency,
    compare_standard_regimens,
    evaluate_therapy_dose_matrix,
    format_markdown_report,
    simulate_efferocytosis_dynamics,
)


class TestMacrophageEfferocytosisModel(unittest.TestCase):
    """Test suite for Alveolar & Airway Macrophage Efferocytosis Impairment & NET Degradation Model."""

    def test_01_trial_provenance_and_snapshot_grounding(self):
        """Clinical trial IDs cited must exist in tools/trial_atlas/snapshot.json."""
        snapshot_file = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
        self.assertTrue(snapshot_file.exists(), "snapshot.json must exist in tools/trial_atlas/")

        with open(snapshot_file, "r", encoding="utf-8") as f:
            snapshot_data = json.load(f)

        studies = snapshot_data.get("studies", [])
        study_ncts = {s["nct_id"] for s in studies if "nct_id" in s}

        # Check rhDNase trial and antioxidant / NAC / GSH trials
        self.assertIn("NCT00843817", study_ncts, "RhDNase NCT00843817 must be in snapshot.json")
        self.assertIn("NCT02088216", study_ncts, "NAC NCT02088216 must be in snapshot.json")
        self.assertIn("NCT00506688", study_ncts, "Glutathione NCT00506688 must be in snapshot.json")

    def test_02_phagolysosomal_acidification_efficiency(self):
        """Test pH-dependent phagolysosomal digestive efficiency eta_phago(pH)."""
        eta_wt = calculate_phago_efficiency(DEFAULT_PH_WT, DEFAULT_PK_LYS)
        eta_cf = calculate_phago_efficiency(DEFAULT_PH_CF, DEFAULT_PK_LYS)

        # At WT pH 4.8 (acidic): eta should be approx 0.715
        self.assertAlmostEqual(eta_wt, 0.715, places=2)

        # At CF pH 6.2 (blunted): eta should be approx 0.091
        self.assertAlmostEqual(eta_cf, 0.091, places=2)

        # Verification of the ~8-fold digestive impairment documented in Di et al. 2006
        ratio = eta_wt / eta_cf
        self.assertGreater(ratio, 7.0)
        self.assertLess(ratio, 9.0)

    def test_03_mertk_scavenger_cleavage_steady_state(self):
        """Test active NE-mediated cleavage of MerTK scavenger receptors."""
        params_no_ne = EfferocytosisParams(ne_conc=0.0)
        mertk_no_ne = calculate_mertk_steady_state(params_no_ne)
        self.assertAlmostEqual(mertk_no_ne, DEFAULT_MERTK_0_RECEPTORS_CELL, places=1)

        # With 180 nM active NE in CF sputum, surface MerTK should be depleted > 90%
        params_cf = EfferocytosisParams(ne_conc=DEFAULT_NE_CONC_NM)
        mertk_cf = calculate_mertk_steady_state(params_cf)
        self.assertLess(mertk_cf, 0.10 * DEFAULT_MERTK_0_RECEPTORS_CELL)

        # NE inhibitor (e.g. 90% inhibition) rescues intact MerTK
        params_inh = EfferocytosisParams(ne_conc=DEFAULT_NE_CONC_NM, ne_inhibitor_potency=0.90)
        mertk_inh = calculate_mertk_steady_state(params_inh)
        self.assertGreater(mertk_inh, mertk_cf * 5.0)

    def test_04_efferocytosis_clearance_kinetics(self):
        """Test instantaneous efferocytosis rate calculation and boundary conditions."""
        params = EfferocytosisParams()
        eta = calculate_phago_efficiency(DEFAULT_PH_WT)
        mertk = DEFAULT_MERTK_0_RECEPTORS_CELL

        # Zero apoptotic cells -> zero clearance
        self.assertEqual(calculate_efferocytosis_rate(0.0, mertk, eta, params), 0.0)

        # Zero intact receptors -> zero clearance
        self.assertEqual(calculate_efferocytosis_rate(1e5, 0.0, eta, params), 0.0)

        # Zero digestive efficiency -> zero clearance
        self.assertEqual(calculate_efferocytosis_rate(1e5, mertk, 0.0, params), 0.0)

        # Normal physiological capacity
        rate = calculate_efferocytosis_rate(DEFAULT_KM_EFF_CELLS_ML, mertk, eta, params)
        expected_rate = params.k_eff_0 * eta * 1.0 * params.m_phi * 0.5
        self.assertAlmostEqual(rate, expected_rate, places=1)

    def test_05_untreated_cf_secondary_necrosis_and_net_accumulation(self):
        """Test that untreated CF dynamics lead to secondary necrotic accumulation and NET hydrogel."""
        params = EfferocytosisParams()
        res = simulate_efferocytosis_dynamics(params, time_hours=48.0)

        # Apoptotic PMNs unengulfed progress to secondary necrosis
        self.assertGreater(res.final_necrotic_pmn, 5e4)
        # Low efferocytosis clearance fraction in CF (< 25%)
        self.assertLess(res.efferocytosis_clearance_fraction, 0.25)
        # Elevated NET eDNA mass (> 5 ug/mL)
        self.assertGreater(res.final_net_mass, 5.0)
        # Toxic extracellular histones present (> 2 ug/mL)
        self.assertGreater(res.final_histone_mass, 2.0)

    def test_06_cftr_and_ne_inhibitor_synergy(self):
        """Test synergy between NE inhibition (receptor protection) and CFTR rescue (acidification)."""
        regimens = compare_standard_regimens(time_hours=48.0)
        cf = regimens["Untreated CF Airway"]
        nei = regimens["NE Inhibitor (85% Block)"]
        cftr = regimens["CFTR Rescue (pH 4.8)"]
        triple = regimens["Triple Synergy Combination"]
        norm = regimens["Normalized Healthy Baseline"]

        # Both monotherapies improve efferocytosis clearance over untreated CF
        self.assertGreater(nei.efferocytosis_clearance_fraction, cf.efferocytosis_clearance_fraction)
        self.assertGreater(cftr.efferocytosis_clearance_fraction, cf.efferocytosis_clearance_fraction)

        # Combination therapy achieves superior clearance fraction over either monotherapy
        self.assertGreater(triple.efferocytosis_clearance_fraction, nei.efferocytosis_clearance_fraction)
        self.assertGreater(triple.efferocytosis_clearance_fraction, cftr.efferocytosis_clearance_fraction)

        # Multi-fold synergy over monotherapies (> 5x higher than either monotherapy)
        self.assertGreater(triple.efferocytosis_clearance_fraction, 5.0 * nei.efferocytosis_clearance_fraction)
        self.assertGreater(triple.efferocytosis_clearance_fraction, 5.0 * cftr.efferocytosis_clearance_fraction)

        # Normalized healthy baseline achieves physiologic efferocytic clearance (> 70%)
        self.assertGreater(norm.efferocytosis_clearance_fraction, 0.70)

        # Combination therapy and healthy baseline reduce toxic histone accumulation
        self.assertLess(triple.final_histone_mass, cf.final_histone_mass)
        self.assertLess(norm.final_histone_mass, cf.final_histone_mass * 0.40)

    def test_07_rhdnase_net_degradation(self):
        """Test dornase alfa (rhDNase I) enhancement of NET eDNA degradation (NCT00843817)."""
        regimens = compare_standard_regimens(time_hours=48.0)
        cf = regimens["Untreated CF Airway"]
        dnase = regimens["rhDNase I (Dornase Alfa)"]

        # Dornase alfa accelerates clearance, dropping final NET mass by > 80%
        self.assertLess(dnase.final_net_mass, cf.final_net_mass * 0.20)

    def test_08_cli_and_disclaimers(self):
        """Test CLI flags and mandatory NLM registry disclaimer presence."""
        script_path = REPO_ROOT / "tools" / "discovery" / "macrophage_efferocytosis_model.py"

        # Test --compare flag
        result = subprocess.run(
            [sys.executable, str(script_path), "--compare", "--hours", "12.0"],
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertIn("ClinicalTrials.gov does not evaluate safety or efficacy", result.stdout)
        self.assertIn("NCT00843817", result.stdout)
        self.assertIn("Triple Synergy Combination", result.stdout)

        # Test --dose-matrix flag
        result_matrix = subprocess.run(
            [sys.executable, str(script_path), "--dose-matrix"],
            capture_output=True,
            text=True,
            check=True,
        )
        matrix_data = json.loads(result_matrix.stdout)
        self.assertIsInstance(matrix_data, list)
        self.assertGreater(len(matrix_data), 10)


if __name__ == "__main__":
    unittest.main()
