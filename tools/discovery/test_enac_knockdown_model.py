#!/usr/bin/env python3
"""Unit tests for enac_knockdown_model.py (T-0200)."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from tools.discovery.enac_knockdown_model import (
    ENAC_DRUG_CATALOG,
    NLM_REGISTRY_DISCLAIMER,
    MEDICAL_DISCLAIMER,
    ASL_HEIGHT_COLLAPSED_UM,
    ASL_HEIGHT_PHYSIOLOGICAL_UM,
    calculate_hill_inhibition,
    calculate_asl_height_um,
    simulate_enac_knockdown,
    compare_all_enac_drugs,
)


class TestENaCKnockdownModel(unittest.TestCase):
    """Test suite for ENaC/SCNN1A RNAi vs small-molecule knockdown and ASL rehydration model."""

    def test_01_catalog_structure_and_nct_provenance(self):
        """All 4 modalities must be defined and clinical trial NCTs verified in snapshot.json."""
        expected_keys = {"ARO-ENaC", "BI-1265162", "Amiloride", "CF-Vehicle"}
        self.assertEqual(set(ENAC_DRUG_CATALOG.keys()), expected_keys)

        snapshot_file = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
        self.assertTrue(snapshot_file.exists())
        with open(snapshot_file, "r", encoding="utf-8") as f:
            snapshot_data = json.load(f)

        studies = snapshot_data.get("studies", [])
        study_ncts = {s["nct_id"] for s in studies if "nct_id" in s}

        self.assertIn("NCT04375514", study_ncts)
        self.assertIn("NCT02343445", study_ncts)
        self.assertEqual(ENAC_DRUG_CATALOG["ARO-ENaC"].nct_id, "NCT04375514")
        self.assertEqual(ENAC_DRUG_CATALOG["BI-1265162"].nct_id, "NCT02343445")

    def test_02_hill_inhibition_boundary_conditions(self):
        """Hill inhibition must obey standard pharmacodynamic limits."""
        # Concentration 0 -> 0% inhibition
        self.assertEqual(calculate_hill_inhibition(0.0, 10.0, 1.0), 0.0)

        # Concentration == IC50 -> 50% inhibition
        self.assertAlmostEqual(calculate_hill_inhibition(10.0, 10.0, 1.0), 0.5, places=5)

        # High concentration -> approaches 100%
        self.assertGreater(calculate_hill_inhibition(10000.0, 10.0, 1.0), 0.99)

        # Monotonicity
        prev = 0.0
        for c in [0.1, 1.0, 5.0, 10.0, 50.0, 200.0]:
            cur = calculate_hill_inhibition(c, 10.0, 1.0)
            self.assertGreater(cur, prev)
            prev = cur

    def test_03_asl_height_rehydration_curve(self):
        """ASL height must map 0% inhibition to collapsed 3.2 um and 100% to >= 7.0 um."""
        j_na_0, h_0 = calculate_asl_height_um(0.0)
        self.assertAlmostEqual(j_na_0, 4.5, places=2)
        self.assertAlmostEqual(h_0, ASL_HEIGHT_COLLAPSED_UM, places=2)

        j_na_1, h_1 = calculate_asl_height_um(1.0)
        self.assertLess(j_na_1, 1.5)
        self.assertGreaterEqual(h_1, ASL_HEIGHT_PHYSIOLOGICAL_UM)

        # Partial inhibition intermediate height
        _, h_half = calculate_asl_height_um(0.5)
        self.assertGreater(h_half, h_0)
        self.assertLess(h_half, h_1)

    def test_04_aro_enac_sirna_pharmacokinetics(self):
        """ARO-ENaC siRNA must demonstrate prolonged lung tissue persistence (t1/2 > 300 hours)."""
        cfg = ENAC_DRUG_CATALOG["ARO-ENaC"]
        self.assertEqual(cfg.target_gene, "SCNN1A")
        self.assertEqual(cfg.modality_class, "rnai_sirna")
        self.assertGreaterEqual(cfg.lung_half_life_hours, 400.0)

        profs, rep = simulate_enac_knockdown("ARO-ENaC", total_hours=504.0)
        self.assertEqual(rep.drug_id, "ARO-ENaC")
        self.assertGreater(rep.average_inhibition_pct, 80.0)
        self.assertGreater(rep.average_asl_height_um, 6.5)

    def test_05_renal_therapeutic_index_contrast(self):
        """ARO-ENaC must have wide renal safety margin (TI >= 50) while small molecules have narrow TI."""
        _, rep_sirna = simulate_enac_knockdown("ARO-ENaC")
        _, rep_bi = simulate_enac_knockdown("BI-1265162")
        _, rep_ami = simulate_enac_knockdown("Amiloride")

        self.assertGreaterEqual(rep_sirna.renal_therapeutic_index, 50.0)
        self.assertEqual(rep_sirna.hyperkalemia_risk_level, "LOW")

        self.assertLessEqual(rep_bi.renal_therapeutic_index, 3.0)
        self.assertLessEqual(rep_ami.renal_therapeutic_index, 2.0)

    def test_06_bi_1265162_small_molecule_profile(self):
        """BI-1265162 must exhibit short airway half-life (~2.2h) and require TID dosing."""
        cfg = ENAC_DRUG_CATALOG["BI-1265162"]
        self.assertEqual(cfg.nct_id, "NCT02343445")
        self.assertEqual(cfg.lung_half_life_hours, 2.2)
        self.assertEqual(cfg.standard_interval_hours, 8.0)

        profs, rep = simulate_enac_knockdown("BI-1265162", total_hours=24.0)
        self.assertGreater(rep.peak_inhibition_pct, 95.0)
        self.assertGreater(rep.peak_asl_height_um, 7.0)

    def test_07_amiloride_benchmark(self):
        """Amiloride benchmark must have very short half-life (<1 hour) and high systemic escape."""
        cfg = ENAC_DRUG_CATALOG["Amiloride"]
        self.assertLessEqual(cfg.lung_half_life_hours, 1.0)
        self.assertGreaterEqual(cfg.systemic_escape_fraction, 0.50)

        profs, rep = simulate_enac_knockdown("Amiloride", total_hours=12.0)
        self.assertEqual(rep.hyperkalemia_risk_level, "HIGH")

    def test_08_cf_vehicle_control_baseline(self):
        """CF-Vehicle must produce 0% inhibition, collapsed ASL (3.2 um), and no renal risk."""
        profs, rep = simulate_enac_knockdown("CF-Vehicle", total_hours=24.0)
        self.assertEqual(rep.peak_inhibition_pct, 0.0)
        self.assertEqual(rep.average_inhibition_pct, 0.0)
        self.assertAlmostEqual(rep.average_asl_height_um, 3.2, places=1)
        self.assertEqual(rep.fraction_time_hydrated, 0.0)

    def test_09_trajectory_non_negativity_and_monotonic_decay(self):
        """Concentrations and ASL height must remain physically bounded throughout simulation."""
        profs, _ = simulate_enac_knockdown("ARO-ENaC", total_hours=100.0)
        for p in profs:
            self.assertGreaterEqual(p.lung_conc_nm, 0.0)
            self.assertGreaterEqual(p.systemic_serum_conc_nm, 0.0)
            self.assertGreaterEqual(p.enac_inhibition_fraction, 0.0)
            self.assertLessEqual(p.enac_inhibition_fraction, 1.0)
            self.assertGreaterEqual(p.asl_height_um, 3.2)
            self.assertLessEqual(p.asl_height_um, 8.0)

    def test_10_simulation_report_consistency(self):
        """Simulation report metrics must maintain peak >= average >= trough relationships."""
        for drug_key in ENAC_DRUG_CATALOG:
            _, rep = simulate_enac_knockdown(drug_key, total_hours=48.0)
            self.assertGreaterEqual(rep.peak_inhibition_pct, rep.average_inhibition_pct)
            self.assertGreaterEqual(rep.average_inhibition_pct, rep.trough_inhibition_pct - 0.1)
            self.assertGreaterEqual(rep.peak_asl_height_um, rep.average_asl_height_um)
            self.assertGreaterEqual(rep.average_asl_height_um, rep.trough_asl_height_um - 0.1)

    def test_11_compare_all_enac_drugs(self):
        """compare_all_enac_drugs must return records for all 4 modalities."""
        comp = compare_all_enac_drugs(total_hours=72.0)
        self.assertEqual(len(comp), 4)
        for key in ENAC_DRUG_CATALOG:
            self.assertIn(key, comp)
            self.assertIn("report", comp[key])
            self.assertIn("sample_profiles", comp[key])

    def test_12_cli_json_and_self_test(self):
        """CLI invocation with --json, --compare-all, and --self-test must succeed."""
        script_path = REPO_ROOT / "tools" / "discovery" / "enac_knockdown_model.py"

        # Check self-test
        res_self = subprocess.run([sys.executable, str(script_path), "--self-test"], capture_output=True, text=True, check=True)
        self.assertIn("ALL ENAC KNOCKDOWN MODEL INVARIANTS PASSED", res_self.stdout)

        # Check json mode
        res_json = subprocess.run([sys.executable, str(script_path), "--drug", "ARO-ENaC", "--json"], capture_output=True, text=True, check=True)
        data = json.loads(res_json.stdout)
        self.assertEqual(data["report"]["nct_id"], "NCT04375514")

        # Check compare-all json
        res_comp = subprocess.run([sys.executable, str(script_path), "--compare-all", "--json"], capture_output=True, text=True, check=True)
        data_comp = json.loads(res_comp.stdout)
        self.assertEqual(len(data_comp), 4)

    def test_13_disclaimers_and_anti_hallucination_guardrails(self):
        """Must include NLM and medical disclaimers and omit prohibited marketing phrases."""
        self.assertIn("attribution to the U.S. National Library of Medicine", NLM_REGISTRY_DISCLAIMER)
        self.assertEqual(MEDICAL_DISCLAIMER, "Research, not medical advice.")

        script_path = REPO_ROOT / "tools" / "discovery" / "enac_knockdown_model.py"
        with open(script_path, "r", encoding="utf-8") as f:
            content = f.read()

        prohibited = [
            "proven to cure",
            "guarantees clinical benefit",
            "recommended treatment regimen",
            "safe and effective for patients",
        ]
        for phrase in prohibited:
            self.assertNotIn(phrase.lower(), content.lower())


if __name__ == "__main__":
    unittest.main()
