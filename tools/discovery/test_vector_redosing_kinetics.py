#!/usr/bin/env python3
"""Unit tests for vector_redosing_kinetics.py (T-0198)."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from tools.discovery.vector_redosing_kinetics import (
    MODALITY_CATALOG,
    NLM_REGISTRY_DISCLAIMER,
    MEDICAL_DISCLAIMER,
    THERAPEUTIC_THRESHOLD_CFTR_PCT,
    calculate_anamnestic_titer,
    calculate_hill_delivery_efficiency,
    simulate_redosing_trajectory,
    compare_all_modalities,
)


class TestVectorRedosingKinetics(unittest.TestCase):
    """Test suite for inhaled vector repeat-dosing and epithelial turnover decay model."""

    def test_01_catalog_structure_and_nct_provenance(self):
        """All 5 vector modalities must be defined and their NCT IDs exist in snapshot.json."""
        expected_modalities = {"4D-710", "BI-3720931", "KB407", "VX-522", "RCT2100"}
        self.assertEqual(set(MODALITY_CATALOG.keys()), expected_modalities)

        snapshot_file = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
        self.assertTrue(snapshot_file.exists(), f"snapshot.json missing at {snapshot_file}")
        with open(snapshot_file, "r", encoding="utf-8") as f:
            snapshot_data = json.load(f)

        studies = snapshot_data.get("studies", [])
        study_ncts = {s["nct_id"] for s in studies if "nct_id" in s}
        for mod_id, cfg in MODALITY_CATALOG.items():
            self.assertTrue(cfg.nct_id.startswith("NCT"))
            self.assertIn(
                cfg.nct_id,
                study_ncts,
                f"Modality {mod_id} NCT ID {cfg.nct_id} not found in trial atlas snapshot",
            )
            self.assertGreater(cfg.luminal_half_life_days, 0.0)
            self.assertGreater(cfg.ic50_titer, 0.0)
            self.assertGreater(cfg.hill_coefficient, 0.0)

    def test_02_hill_neutralization_limits(self):
        """Hill delivery efficiency must obey mathematical boundary conditions."""
        # Titer 0 -> efficiency 100%, neutralization 0%
        eff_0, neut_0 = calculate_hill_delivery_efficiency(0.0, 1.0, 50.0, 1.5)
        self.assertAlmostEqual(eff_0, 1.0)
        self.assertAlmostEqual(neut_0, 0.0)

        # Titer == IC50 -> efficiency 50%, neutralization 50%
        eff_mid, neut_mid = calculate_hill_delivery_efficiency(50.0, 1.0, 50.0, 1.5)
        self.assertAlmostEqual(eff_mid, 0.5)
        self.assertAlmostEqual(neut_mid, 0.5)

        # Titer >> IC50 -> efficiency -> 0%
        eff_high, neut_high = calculate_hill_delivery_efficiency(5000.0, 1.0, 50.0, 1.5)
        self.assertLess(eff_high, 0.001)
        self.assertGreater(neut_high, 0.999)

        # Monotonic decrease with rising titer
        prev_eff = 1.0
        for t in [1.0, 10.0, 50.0, 100.0, 500.0]:
            cur_eff, _ = calculate_hill_delivery_efficiency(t, 1.0, 50.0, 1.5)
            self.assertLess(cur_eff, prev_eff)
            prev_eff = cur_eff

    def test_03_anamnestic_titer_accumulation(self):
        """Titer should accumulate anamnestically with doses and decay exponentially."""
        dose_times = [0.0, 30.0, 60.0]
        # At day 0 (pre-dose 1): baseline titer
        t0 = calculate_anamnestic_titer(dose_times[:0], 0.0, baseline_titer=2.0, boost_per_dose=100.0, decay_half_life_days=60.0)
        self.assertEqual(t0, 2.0)

        # At day 30 (pre-dose 2): baseline + decayed boost from dose 1
        t30 = calculate_anamnestic_titer(dose_times[:1], 30.0, baseline_titer=2.0, boost_per_dose=100.0, decay_half_life_days=60.0)
        expected_boost_1 = 100.0 * (0.5 ** (30.0 / 60.0))
        self.assertAlmostEqual(t30, 2.0 + expected_boost_1, places=3)

        # At day 60 (pre-dose 3): baseline + decayed dose 1 + decayed dose 2
        t60 = calculate_anamnestic_titer(dose_times[:2], 60.0, baseline_titer=2.0, boost_per_dose=100.0, decay_half_life_days=60.0)
        self.assertGreater(t60, t30)

    def test_04_aav_vs_mrna_attenuation_contrast(self):
        """AAV 4D-710 must show severe redosing attenuation, whereas mRNA VX-522 maintains >85%."""
        _, _, sum_aav = simulate_redosing_trajectory("4D-710", dose_count=4, interval_days=30.0)
        _, _, sum_mrna = simulate_redosing_trajectory("VX-522", dose_count=4, interval_days=30.0)

        self.assertLess(sum_aav.redosing_attenuation_ratio, 0.20)
        self.assertGreater(sum_mrna.redosing_attenuation_ratio, 0.85)
        self.assertGreater(sum_mrna.final_dose_efficiency, sum_aav.final_dose_efficiency)

    def test_05_integrative_lentivirus_long_term_durability(self):
        """BI-3720931 lentivirus must maintain higher expression at Day 365 due to basal replenishment."""
        _, traj_lenti, sum_lenti = simulate_redosing_trajectory("BI-3720931", dose_count=2, interval_days=30.0, total_days=365.0)
        _, traj_aav, sum_aav = simulate_redosing_trajectory("4D-710", dose_count=2, interval_days=30.0, total_days=365.0)

        self.assertTrue(sum_lenti.is_curative_integrative)
        self.assertFalse(sum_aav.is_curative_integrative)
        final_lenti = traj_lenti[-1].total_cftr_pct_wt
        final_aav = traj_aav[-1].total_cftr_pct_wt
        self.assertGreater(final_lenti, final_aav)
        self.assertGreater(traj_lenti[-1].basal_stem_pct, traj_aav[-1].basal_stem_pct)

    def test_06_kb407_hsv1_immunology(self):
        """KB407 must reflect high baseline NAb titer and strong boost on redosing."""
        cfg = MODALITY_CATALOG["KB407"]
        self.assertEqual(cfg.baseline_nab_titer, 5.0)
        self.assertGreater(cfg.nab_boost_per_dose, 300.0)

        doses, _, summary = simulate_redosing_trajectory("KB407", dose_count=3, interval_days=28.0)
        self.assertEqual(len(doses), 3)
        self.assertLess(summary.redosing_attenuation_ratio, 0.35)

    def test_07_rct2100_mrna_repeat_dosing(self):
        """RCT2100 must support repeat dosing with high efficiency retention."""
        cfg = MODALITY_CATALOG["RCT2100"]
        self.assertEqual(cfg.genome_type, "transient_mrna")
        self.assertEqual(cfg.nct_id, "NCT06237335")

        doses, traj, summary = simulate_redosing_trajectory("RCT2100", dose_count=6, interval_days=14.0)
        self.assertGreater(summary.redosing_attenuation_ratio, 0.80)
        self.assertGreater(summary.peak_cftr_pct_wt, 30.0)

    def test_08_trajectory_continuity_and_non_negativity(self):
        """All simulated CFTR percentages and titers must remain strictly non-negative."""
        for mod_key in MODALITY_CATALOG:
            _, traj, summary = simulate_redosing_trajectory(mod_key, dose_count=3, interval_days=30.0)
            for pt in traj:
                self.assertGreaterEqual(pt.ciliated_luminal_pct, 0.0)
                self.assertGreaterEqual(pt.basal_stem_pct, 0.0)
                self.assertGreaterEqual(pt.total_cftr_pct_wt, 0.0)
                self.assertGreaterEqual(pt.active_nab_titer, 0.0)
                self.assertAlmostEqual(pt.total_cftr_pct_wt, pt.ciliated_luminal_pct + pt.basal_stem_pct, places=2)

            self.assertGreaterEqual(summary.peak_cftr_pct_wt, summary.trough_cftr_pct_wt)
            self.assertGreaterEqual(summary.time_above_15pct_threshold_days, 0.0)

    def test_09_compare_all_modalities_coverage(self):
        """compare_all_modalities must return complete records for all 5 modalities."""
        comp = compare_all_modalities(dose_count=3, interval_days=30.0, total_days=120.0)
        self.assertEqual(len(comp), 5)
        for key in MODALITY_CATALOG:
            self.assertIn(key, comp)
            self.assertIn("summary", comp[key])
            self.assertIn("doses", comp[key])
            self.assertIn("sample_trajectory_points", comp[key])
            self.assertEqual(comp[key]["summary"]["nct_id"], MODALITY_CATALOG[key].nct_id)

    def test_10_input_validation(self):
        """Invalid inputs must raise ValueError or KeyError."""
        with self.assertRaises(KeyError):
            simulate_redosing_trajectory("NON_EXISTENT_VECTOR")
        with self.assertRaises(ValueError):
            simulate_redosing_trajectory("4D-710", dose_count=0)
        with self.assertRaises(ValueError):
            simulate_redosing_trajectory("4D-710", interval_days=0.5)

    def test_11_cli_json_and_compare_execution(self):
        """Executing CLI with --json and --compare-all must return exit code 0 and valid JSON."""
        script_path = REPO_ROOT / "tools" / "discovery" / "vector_redosing_kinetics.py"
        res_json = subprocess.run(
            [sys.executable, str(script_path), "--modality", "4D-710", "--doses", "2", "--json"],
            capture_output=True,
            text=True,
            check=True,
        )
        data = json.loads(res_json.stdout)
        self.assertIn("summary", data)
        self.assertEqual(data["summary"]["nct_id"], "NCT05248230")

        res_comp = subprocess.run(
            [sys.executable, str(script_path), "--compare-all", "--json"],
            capture_output=True,
            text=True,
            check=True,
        )
        data_comp = json.loads(res_comp.stdout)
        self.assertEqual(len(data_comp), 5)

    def test_12_cli_self_test(self):
        """CLI --self-test flag must exit with code 0."""
        script_path = REPO_ROOT / "tools" / "discovery" / "vector_redosing_kinetics.py"
        res = subprocess.run(
            [sys.executable, str(script_path), "--self-test"],
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertIn("ALL VECTOR REDOSING KINETICS INVARIANTS PASSED", res.stdout)

    def test_13_disclaimers_and_anti_hallucination_guardrails(self):
        """Must include NLM and medical disclaimers and omit prohibited overclaim language."""
        self.assertIn("attribution to the U.S. National Library of Medicine", NLM_REGISTRY_DISCLAIMER)
        self.assertEqual(MEDICAL_DISCLAIMER, "Research, not medical advice.")

        with open(REPO_ROOT / "tools" / "discovery" / "vector_redosing_kinetics.py", "r", encoding="utf-8") as f:
            content = f.read()

        prohibited = [
            "proven to cure",
            "guarantees clinical benefit",
            "recommended treatment regimen",
            "safe and effective for patients",
        ]
        for phrase in prohibited:
            self.assertNotIn(phrase.lower(), content.lower(), f"Prohibited phrase '{phrase}' found in script")


if __name__ == "__main__":
    unittest.main()
