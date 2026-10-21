#!/usr/bin/env python3
"""Unit tests for phage_omko1_tradeoff.py (T-0246)."""

from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from tools.discovery.phage_omko1_tradeoff import (
    DEFAULT_B_R_0_CFU_ML,
    DEFAULT_B_S_0_CFU_ML,
    DEFAULT_BURST_SIZE,
    DEFAULT_C_FIT,
    DEFAULT_CIPRO_CONC_UG_ML,
    DEFAULT_DELTA_P_PER_MIN,
    DEFAULT_E_MAX_PER_MIN,
    DEFAULT_HILL_COEFF,
    DEFAULT_K_ADS_ML_PER_CFU_MIN,
    DEFAULT_K_CFU_ML,
    DEFAULT_LATENT_PERIOD_MIN,
    DEFAULT_MIC_R_UG_ML,
    DEFAULT_MIC_S_UG_ML,
    DEFAULT_PHAGE_0_PFU_ML,
    DEFAULT_R_S_PER_MIN,
    NLM_REGISTRY_DISCLAIMER,
    RESEARCH_BOUNDARY,
    PhageTradeoffParams,
    calculate_hill_kill_rate,
    compare_standard_regimens,
    evaluate_synergy_dose_matrix,
    format_markdown_report,
    simulate_omko1_tradeoff_dynamics,
)


class TestPhageOMKO1Tradeoff(unittest.TestCase):
    """Test suite for Phage OMKO1 Evolutionary Trade-off and Resensitization Dynamics."""

    def test_01_trial_provenance_and_snapshot_grounding(self):
        """Clinical trial IDs cited must exist in tools/trial_atlas/snapshot.json."""
        snapshot_file = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
        self.assertTrue(snapshot_file.exists(), "snapshot.json must exist in tools/trial_atlas/")

        with open(snapshot_file, "r", encoding="utf-8") as f:
            snapshot_data = json.load(f)

        studies = snapshot_data.get("studies", [])
        study_ncts = {s["nct_id"] for s in studies if "nct_id" in s}

        # Check CYPHY Yale trial and Inhaled Phage trial
        self.assertIn("NCT04684641", study_ncts, "CYPHY NCT04684641 must be in snapshot.json")
        self.assertIn("NCT05010577", study_ncts, "Inhaled Phage NCT05010577 must be in snapshot.json")

    def test_02_hill_kill_rate_and_resensitization_ratio(self):
        """Pharmacodynamic Hill kill rate and MIC fold drop must adhere to biophysical models."""
        params = PhageTradeoffParams()
        # Fold reduction: 32.0 / 1.0 = 32.0x (must be >= 16x)
        self.assertGreaterEqual(params.resensitization_fold_drop, 16.0)
        self.assertAlmostEqual(params.resensitization_fold_drop, 32.0, places=3)

        # At cipro = 10 ug/mL:
        # For B_S (MIC = 32 ug/mL), kill rate should be minimal
        kill_s = calculate_hill_kill_rate(10.0, 32.0, DEFAULT_E_MAX_PER_MIN, DEFAULT_HILL_COEFF)
        # For B_R (MIC = 1.0 ug/mL), kill rate should approach E_max
        kill_r = calculate_hill_kill_rate(10.0, 1.0, DEFAULT_E_MAX_PER_MIN, DEFAULT_HILL_COEFF)

        self.assertLess(kill_s, 0.015)
        self.assertGreater(kill_r, 0.075)
        self.assertGreater(kill_r / kill_s, 5.0)

        # Zero concentration produces zero kill
        self.assertEqual(calculate_hill_kill_rate(0.0, 1.0), 0.0)

    def test_03_parameter_validation_and_bounds(self):
        """Parameter container must validate physiological boundaries."""
        # Valid parameters should pass
        valid_p = PhageTradeoffParams()
        valid_p.validate()

        # Negative growth rate must raise ValueError
        with self.assertRaises(ValueError):
            PhageTradeoffParams(r_s=-0.01).validate()

        # Fitness cost >= 1.0 must raise ValueError
        with self.assertRaises(ValueError):
            PhageTradeoffParams(c_fit=1.2).validate()

        # Inverted MICs (MIC_S <= MIC_R) must raise ValueError
        with self.assertRaises(ValueError):
            PhageTradeoffParams(mic_s=1.0, mic_r=32.0).validate()

    def test_04_untreated_control_logistic_growth(self):
        """Untreated control must exhibit logistic growth toward carrying capacity without eradication."""
        params = PhageTradeoffParams(cipro_conc_ug_ml=0.0, phage_0=0.0)
        res = simulate_omko1_tradeoff_dynamics(params, duration_hours=12.0, scenario_name="Untreated")

        self.assertFalse(res.eradication_achieved)
        self.assertGreater(res.final_b_tot, res.initial_b_tot)
        self.assertLessEqual(res.final_b_tot, params.carrying_capacity * 1.05)
        self.assertEqual(res.final_phage, 0.0)

    def test_05_monotherapy_failures_and_resistance_escape(self):
        """Monotherapies must fail to sterilize: Cipro fails on B_S; Phage selects for B_R escape."""
        regimens = compare_standard_regimens(duration_hours=24.0)

        cipro_res = regimens["cipro_monotherapy"]
        phage_res = regimens["phage_monotherapy"]

        # Cipro monotherapy fails to eradicate due to high MIC_S
        self.assertFalse(cipro_res.eradication_achieved)
        self.assertGreater(cipro_res.final_b_tot, 1e6)

        # Phage monotherapy drives high selection for resistant mutant (B_R fraction > 90%)
        self.assertGreater(phage_res.final_fraction_resistant, 0.85)
        self.assertGreater(phage_res.final_b_r, DEFAULT_B_R_0_CFU_ML)

    def test_06_dual_evolutionary_trap_synergy(self):
        """Dual therapy (OMKO1 + ciprofloxacin) must trigger synergistic collapse of both strains."""
        regimens = compare_standard_regimens(duration_hours=24.0)
        dual_res = regimens["dual_evolutionary_trap"]

        # Dual therapy achieves dramatic log reduction (> 5 logs) or complete sterile eradication
        self.assertTrue(dual_res.eradication_achieved or dual_res.log10_b_tot_reduction >= 5.0)
        self.assertGreater(dual_res.evolutionary_vice_score, 4.0)

        # Vice score of dual therapy must exceed both monotherapies
        self.assertGreater(dual_res.evolutionary_vice_score, regimens["cipro_monotherapy"].evolutionary_vice_score)
        self.assertGreater(dual_res.evolutionary_vice_score, regimens["phage_monotherapy"].evolutionary_vice_score)

    def test_07_dose_matrix_evaluation(self):
        """Dose matrix grid sweep must execute across defined ranges."""
        matrix = evaluate_synergy_dose_matrix(
            cipro_levels=[0.0, 5.0, 10.0],
            phage_levels=[0.0, 1e7, 1e8],
            duration_hours=12.0,
        )
        self.assertEqual(matrix["total_simulations"], 9)
        self.assertEqual(len(matrix["entries"]), 9)

        # High cipro + high phage entry must have higher log reduction than 0 + 0
        entry_zero = next(e for e in matrix["entries"] if e["cipro_conc_ug_ml"] == 0.0 and e["phage_dose_pfu_ml"] == 0.0)
        entry_high = next(e for e in matrix["entries"] if e["cipro_conc_ug_ml"] == 10.0 and e["phage_dose_pfu_ml"] == 1e8)

        self.assertGreater(entry_high["log10_reduction"], entry_zero["log10_reduction"])

    def test_08_cli_execution_and_disclaimers(self):
        """CLI execution must output JSON and contain mandatory regulatory disclaimers."""
        script_path = REPO_ROOT / "tools" / "discovery" / "phage_omko1_tradeoff.py"
        proc = subprocess.run(
            [sys.executable, str(script_path), "--duration-hours", "6.0", "--json"],
            capture_output=True,
            text=True,
            check=True,
        )

        data = json.loads(proc.stdout)
        self.assertIn("disclaimer", data)
        self.assertIn("boundary", data)
        self.assertIn("regimens", data)
        self.assertIn("dual_evolutionary_trap", data["regimens"])
        self.assertIn("ClinicalTrials.gov does not evaluate safety or efficacy", data["disclaimer"])


if __name__ == "__main__":
    unittest.main()
