#!/usr/bin/env python3
"""Unit tests for nbd1_thermodynamics.py (T-0199)."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from tools.discovery.nbd1_thermodynamics import (
    PDB_RECORDS,
    VARIANT_CATALOG,
    REGIMEN_CATALOG,
    NLM_REGISTRY_DISCLAIMER,
    MEDICAL_DISCLAIMER,
    calculate_gibbs_helmholtz_unfolding_free_energy,
    evaluate_nbd1_thermodynamics,
    compare_variant_panel,
)


class TestNBD1Thermodynamics(unittest.TestCase):
    """Test suite for CFTR NBD1 thermal unfolding & domain-interface free energy model."""

    def test_01_pdb_and_variant_catalog_integrity(self):
        """PDB cryo-EM structures and all 6 representative CFTR alleles must be defined."""
        self.assertIn("5UAK", PDB_RECORDS)
        self.assertIn("6MSM", PDB_RECORDS)
        self.assertEqual(PDB_RECORDS["5UAK"]["conformation"], "inward-facing, ATP-free")
        self.assertEqual(PDB_RECORDS["6MSM"]["conformation"], "outward-facing, ATP-bound dimer")

        expected_variants = {"WT", "F508del", "R560T", "N1303K", "G85E", "G551D"}
        self.assertEqual(set(VARIANT_CATALOG.keys()), expected_variants)

    def test_02_clinical_trial_nct_provenance_in_snapshot(self):
        """Regimens citing clinical trials must have their NCT IDs verified in snapshot.json."""
        snapshot_file = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
        self.assertTrue(snapshot_file.exists())
        with open(snapshot_file, "r", encoding="utf-8") as f:
            snapshot_data = json.load(f)

        studies = snapshot_data.get("studies", [])
        study_ncts = {s["nct_id"] for s in studies if "nct_id" in s}

        for reg_id, reg_cfg in REGIMEN_CATALOG.items():
            if reg_cfg.nct_id:
                self.assertTrue(reg_cfg.nct_id.startswith("NCT"))
                self.assertIn(
                    reg_cfg.nct_id,
                    study_ncts,
                    f"Regimen {reg_id} NCT ID {reg_cfg.nct_id} missing from snapshot.json",
                )

    def test_03_gibbs_helmholtz_thermodynamics_boundary_conditions(self):
        """Gibbs-Helmholtz unfolding free energy must satisfy physical boundary conditions."""
        tm_c = 49.5
        tm_k = tm_c + 273.15

        # At T = Tm, Delta G_unfold must equal 0.0 (50% unfolded equilibrium)
        dg_tm = calculate_gibbs_helmholtz_unfolding_free_energy(tm_k, tm_k)
        self.assertAlmostEqual(dg_tm, 0.0, places=5)

        # At T < Tm (e.g. 37 deg C), native state is favored (Delta G_unfold > 0)
        t_body_k = 37.0 + 273.15
        dg_body = calculate_gibbs_helmholtz_unfolding_free_energy(t_body_k, tm_k)
        self.assertGreater(dg_body, 0.0)

        # At T > Tm (e.g. 60 deg C), denatured state is favored (Delta G_unfold < 0)
        t_hot_k = 60.0 + 273.15
        dg_hot = calculate_gibbs_helmholtz_unfolding_free_energy(t_hot_k, tm_k)
        self.assertLess(dg_hot, 0.0)

    def test_04_wt_baseline_parameters(self):
        """Wild type CFTR must have Tm = 49.5 deg C, Delta Tm = 0.0, Band C ~ 85%, Po ~ 0.35."""
        wt = evaluate_nbd1_thermodynamics("WT", "None")
        self.assertEqual(wt.effective_tm_celsius, 49.5)
        self.assertEqual(wt.delta_tm_vs_wt_celsius, 0.0)
        self.assertAlmostEqual(wt.band_c_maturation_pct, 85.0, places=1)
        self.assertAlmostEqual(wt.open_probability, 0.35, places=2)
        self.assertAlmostEqual(wt.functional_cftr_pct_wt, 100.0, places=1)

    def test_05_f508del_thermal_destabilization(self):
        """F508del must exhibit severe Delta Tm approx -7.5 deg C, Band C < 2%, and low function."""
        f508 = evaluate_nbd1_thermodynamics("F508del", "None")
        self.assertEqual(f508.effective_tm_celsius, 42.0)
        self.assertEqual(f508.delta_tm_vs_wt_celsius, -7.5)
        self.assertLess(f508.band_c_maturation_pct, 2.0)
        self.assertLess(f508.functional_cftr_pct_wt, 2.0)

    def test_06_r560t_signature_region_instability(self):
        """R560T must exhibit the lowest thermal stability among tested variants (Tm = 40.5 deg C)."""
        r560 = evaluate_nbd1_thermodynamics("R560T", "None")
        f508 = evaluate_nbd1_thermodynamics("F508del", "None")
        self.assertEqual(r560.effective_tm_celsius, 40.5)
        self.assertEqual(r560.delta_tm_vs_wt_celsius, -9.0)
        self.assertLess(r560.effective_tm_celsius, f508.effective_tm_celsius)
        self.assertLess(r560.band_c_maturation_pct, 2.0)

    def test_07_n1303k_and_g85e_domain_distribution(self):
        """N1303K (NBD2) and G85E (TMD1) must have preserved NBD1 Tm despite low maturation."""
        n1303k = evaluate_nbd1_thermodynamics("N1303K", "None")
        g85e = evaluate_nbd1_thermodynamics("G85E", "None")

        self.assertGreaterEqual(n1303k.effective_tm_celsius, 48.5)
        self.assertGreaterEqual(g85e.effective_tm_celsius, 48.0)
        self.assertLess(n1303k.band_c_maturation_pct, 15.0)
        self.assertLess(g85e.band_c_maturation_pct, 5.0)

    def test_08_g551d_gating_defect_and_potentiator_rescue(self):
        """G551D must have normal Tm, normal Band C, low basal Po, and high Ivacaftor rescue."""
        apo = evaluate_nbd1_thermodynamics("G551D", "None")
        self.assertAlmostEqual(apo.effective_tm_celsius, 49.2, places=1)
        self.assertAlmostEqual(apo.band_c_maturation_pct, 80.0, places=1)
        self.assertLess(apo.open_probability, 0.05)
        self.assertLess(apo.functional_cftr_pct_wt, 10.0)

        rescued = evaluate_nbd1_thermodynamics("G551D", "IvacaftorOnly")
        self.assertGreater(rescued.open_probability, 0.25)
        self.assertGreater(rescued.functional_cftr_pct_wt, 70.0)

    def test_09_type1_vs_type3_corrector_stabilization(self):
        """Type III (Elexacaftor) must provide greater direct NBD1 thermal stabilization than Type I (Tezacaftor)."""
        teza = evaluate_nbd1_thermodynamics("F508del", "Tezacaftor")
        elexa = evaluate_nbd1_thermodynamics("F508del", "Elexacaftor")

        self.assertGreater(elexa.effective_tm_celsius, teza.effective_tm_celsius)
        self.assertLess(elexa.folding_delta_g_kcal, teza.folding_delta_g_kcal)

    def test_10_trikafta_cooperative_rescue_nct04043806(self):
        """Trikafta (ETI / NCT04043806) must cooperatively stabilize F508del above 47.0 deg C."""
        trikafta = evaluate_nbd1_thermodynamics("F508del", "Trikafta")
        apo = evaluate_nbd1_thermodynamics("F508del", "None")

        self.assertEqual(trikafta.nct_id, "NCT04043806")
        self.assertGreater(trikafta.effective_tm_celsius, 47.0)
        self.assertGreater(trikafta.band_c_maturation_pct, 60.0)
        self.assertGreater(trikafta.functional_cftr_pct_wt, 40.0)
        self.assertLess(trikafta.folding_delta_g_kcal, apo.folding_delta_g_kcal)

    def test_11_vanza_triple_next_gen_potency_nct06468527(self):
        """VanzaTriple (NCT06468527) must achieve higher Tm boost and favorable Delta G than Trikafta."""
        vanza = evaluate_nbd1_thermodynamics("F508del", "VanzaTriple")
        trikafta = evaluate_nbd1_thermodynamics("F508del", "Trikafta")

        self.assertEqual(vanza.nct_id, "NCT06468527")
        self.assertGreater(vanza.effective_tm_celsius, trikafta.effective_tm_celsius)
        self.assertLess(vanza.folding_delta_g_kcal, trikafta.folding_delta_g_kcal)

    def test_12_compare_variant_panel_coverage(self):
        """compare_variant_panel must return data for all 6 variants with PDB grounding."""
        panel = compare_variant_panel(regimen_key="Trikafta")
        self.assertEqual(len(panel), 6)
        for var_k, data in panel.items():
            self.assertIn("effective_tm_celsius", data)
            self.assertIn("functional_cftr_pct_wt", data)
            self.assertIn("PDB:5UAK", data["pdb_grounding"])
            self.assertIn("PDB:6MSM", data["pdb_grounding"])

    def test_13_cli_execution_and_disclaimers(self):
        """CLI invocation with --json, --compare-variants, and --self-test must succeed without overclaims."""
        script_path = REPO_ROOT / "tools" / "discovery" / "nbd1_thermodynamics.py"

        # Check self-test
        res_self = subprocess.run([sys.executable, str(script_path), "--self-test"], capture_output=True, text=True, check=True)
        self.assertIn("ALL NBD1 THERMODYNAMICS INVARIANTS PASSED", res_self.stdout)

        # Check json mode
        res_json = subprocess.run([sys.executable, str(script_path), "--variant", "F508del", "--json"], capture_output=True, text=True, check=True)
        data = json.loads(res_json.stdout)
        self.assertEqual(data["variant_id"], "F508del")
        self.assertIn("PDB:5UAK", data["pdb_grounding"])

        # Check disclaimers and prohibition of hallucinated marketing claims
        self.assertIn("attribution to the U.S. National Library of Medicine", NLM_REGISTRY_DISCLAIMER)
        self.assertEqual(MEDICAL_DISCLAIMER, "Research, not medical advice.")

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
