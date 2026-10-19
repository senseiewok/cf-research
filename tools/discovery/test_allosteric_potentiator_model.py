#!/usr/bin/env python3
"""Unit tests for allosteric_potentiator_model.py (T-0221)."""
from __future__ import annotations

import json
import math
import subprocess
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from tools.discovery.allosteric_potentiator_model import (
    NLM_REGISTRY_DISCLAIMER,
    MEDICAL_DISCLAIMER,
    PDB_RECORDS,
    LITERATURE_RECORDS,
    MODULATOR_CATALOG,
    REGIMEN_CATALOG,
    VARIANT_BASELINES,
    AllostericModulatorCycle,
    get_pairwise_cooperativity,
    verify_detailed_balance,
    compute_dose_response_surface,
)


class TestAllostericPotentiatorModel(unittest.TestCase):
    """Test suite for CFTR potentiator + triple-corrector cooperative allosteric gating thermodynamic cycle."""

    def test_01_pdb_and_literature_grounding_integrity(self):
        """PDB cryo-EM structures and primary Science literature citation must be defined."""
        self.assertIn("8EIQ", PDB_RECORDS)
        self.assertIn("6MSM", PDB_RECORDS)
        self.assertEqual(PDB_RECORDS["8EIQ"]["pdb_id"], "8EIQ")
        self.assertIn("Ivacaftor", PDB_RECORDS["8EIQ"]["ligands"])
        self.assertIn("Tezacaftor", PDB_RECORDS["8EIQ"]["ligands"])
        self.assertIn("Elexacaftor", PDB_RECORDS["8EIQ"]["ligands"])

        self.assertIn("10.1126/science.ade2216", LITERATURE_RECORDS)
        self.assertIn("Fiedorczuk", LITERATURE_RECORDS["10.1126/science.ade2216"]["citation"])

    def test_02_clinical_trial_nct_provenance_in_snapshot(self):
        """Clinical trials cited in modulator and regimen catalogs must exist in snapshot.json."""
        snapshot_file = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
        self.assertTrue(snapshot_file.exists(), "tools/trial_atlas/snapshot.json must exist")
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

        for mod_id, mod_cfg in MODULATOR_CATALOG.items():
            if mod_cfg.nct_id:
                self.assertTrue(mod_cfg.nct_id.startswith("NCT"))
                self.assertIn(
                    mod_cfg.nct_id,
                    study_ncts,
                    f"Modulator {mod_id} NCT ID {mod_cfg.nct_id} missing from snapshot.json",
                )

    def test_03_modulator_parameter_alpha_efficacy(self):
        """All potentiators and co-potentiators must exhibit alpha = Kd_C / Kd_O > 1.0 (open-state stabilizing)."""
        expected_compounds = {"VX-770", "VX-561", "VX-661", "VX-445", "VX-121"}
        self.assertEqual(set(MODULATOR_CATALOG.keys()), expected_compounds)

        for comp_id, param in MODULATOR_CATALOG.items():
            self.assertGreater(param.alpha, 1.0, f"{comp_id} alpha must be > 1.0")
            self.assertLess(param.ddg_open_kcal, 0.0, f"{comp_id} DeltaDelta G_open must be negative")

        # Ivacaftor (VX-770) and Deutivacaftor (VX-561) are strong potentiators (alpha = 10)
        self.assertAlmostEqual(MODULATOR_CATALOG["VX-770"].alpha, 10.0, places=3)
        self.assertAlmostEqual(MODULATOR_CATALOG["VX-561"].alpha, 10.0, places=3)

        # Elexacaftor (VX-445) has dual activity with alpha ~ 3.0
        self.assertAlmostEqual(MODULATOR_CATALOG["VX-445"].alpha, 3.0, places=3)

        # Vanzacaftor (VX-121) has next-gen dual activity with alpha ~ 3.6
        self.assertAlmostEqual(MODULATOR_CATALOG["VX-121"].alpha, 3.6, places=3)

    def test_04_thermodynamic_detailed_balance_closed_loop(self):
        """Detailed balance around 4-state thermodynamic cycles must hold (cyclic product == 1.0, loop Delta G == 0)."""
        for comp_id in MODULATOR_CATALOG:
            res = verify_detailed_balance(comp_id)
            self.assertTrue(res["detailed_balance_satisfied"], f"{comp_id} detailed balance violated")
            self.assertAlmostEqual(res["cyclic_product"], 1.0, places=5)
            self.assertAlmostEqual(res["loop_free_energy_kcal"], 0.0, places=5)

    def test_05_pairwise_allosteric_cooperativity(self):
        """Pairwise coupling between distinct pockets must exhibit positive cooperativity (gamma > 1.0)."""
        # Type I (VX-661) + Type III (VX-445)
        gamma_661_445 = get_pairwise_cooperativity("VX-661", "VX-445")
        self.assertGreater(gamma_661_445, 1.5, "Tezacaftor + Elexacaftor must exhibit gamma > 1.5")

        # Type III (VX-445) + Potentiator (VX-770)
        gamma_445_770 = get_pairwise_cooperativity("VX-445", "VX-770")
        self.assertGreater(gamma_445_770, 2.0, "Elexacaftor + Ivacaftor hinge coupling must exhibit gamma > 2.0")

        # Next-Gen Type III (VX-121) + Deutivacaftor (VX-561)
        gamma_121_561 = get_pairwise_cooperativity("VX-121", "VX-561")
        self.assertGreater(gamma_121_561, 2.2, "Vanzacaftor + Deutivacaftor coupling must exhibit gamma > 2.2")

    def test_06_trikafta_thermodynamic_synergy_on_f508del(self):
        """Trikafta combination on F508del must produce positive thermodynamic synergy S_123 > 1.0."""
        cycle = AllostericModulatorCycle(variant_id="F508del")
        res = cycle.run_cycle(["VX-445", "VX-661", "VX-770"])

        self.assertGreater(res.p_open_treated, res.p_open_basal)
        self.assertGreater(res.thermodynamic_synergy_s, 1.0, "Trikafta must show positive thermodynamic synergy (S > 1.0)")
        self.assertGreater(res.fold_stimulation, 3.0, "Trikafta should stimulate P_open by > 3-fold on F508del")

    def test_07_nextgen_triple_efficacy_comparison(self):
        """NextGenTriple (VX-121/VX-661/VX-561) must match or exceed Trikafta in P_open and sweat chloride reduction."""
        cycle = AllostericModulatorCycle(variant_id="F508del")
        res_trikafta = cycle.run_cycle(["VX-445", "VX-661", "VX-770"])
        res_nextgen = cycle.run_cycle(["VX-121", "VX-661", "VX-561"])

        self.assertGreaterEqual(res_nextgen.p_open_treated, res_trikafta.p_open_treated * 0.95)
        self.assertLessEqual(res_nextgen.predicted_delta_swcl_mmol_L, res_trikafta.predicted_delta_swcl_mmol_L)
        self.assertLess(res_nextgen.predicted_delta_swcl_mmol_L, -40.0)

    def test_08_single_channel_current_and_driving_force(self):
        """Single-channel current must obey Ohm's law and scale linearly with driving force."""
        cycle = AllostericModulatorCycle(variant_id="WT", membrane_voltage_mV=-60.0, chloride_reversal_mV=0.0)
        p_open = 0.50
        # Conductance 8.4 pS, Vm - ECl = -60 mV = -0.06 V
        # I = 8.4e-12 S * 0.50 * -0.06 V = -0.252 pA
        expected_current_pA = 8.4 * 0.50 * (-0.06)
        curr = cycle.calculate_single_channel_current_pA(p_open)
        self.assertAlmostEqual(curr, expected_current_pA, places=3)
        self.assertLess(curr, 0.0, "Inward chloride current (outward ion flow) is negative at -60 mV")

    def test_09_genotype_response_discrimination(self):
        """G551D must respond robustly to Ivacaftor monotherapy; F508del requires triple combination for high flux."""
        cycle_g551d = AllostericModulatorCycle(variant_id="G551D")
        cycle_f508del = AllostericModulatorCycle(variant_id="F508del")

        # Ivacaftor monotherapy on G551D
        res_g551d_kalydeco = cycle_g551d.run_cycle(["VX-770"])
        self.assertGreater(res_g551d_kalydeco.fold_stimulation, 5.0)

        # Ivacaftor monotherapy on F508del has minimal sweat chloride reduction due to severe lack of surface channels
        res_f508_kalydeco = cycle_f508del.run_cycle(["VX-770"])
        self.assertGreater(res_f508_kalydeco.predicted_delta_swcl_mmol_L, -10.0)

        # Trikafta on F508del produces profound sweat chloride reduction
        res_f508_trikafta = cycle_f508del.run_cycle(["VX-445", "VX-661", "VX-770"])
        self.assertLess(res_f508_trikafta.predicted_delta_swcl_mmol_L, -40.0)

    def test_10_dose_response_monotonicity_and_saturation(self):
        """Increasing potentiator concentration must monotonically increase P_open within [0, 1]."""
        cycle = AllostericModulatorCycle(variant_id="WT")
        concs = [0.001, 0.01, 0.05, 0.20, 1.0, 5.0, 20.0]
        p_opens = [cycle.evaluate_open_probability(["VX-770"], {"VX-770": c}) for c in concs]

        for i in range(len(p_opens) - 1):
            self.assertLessEqual(p_opens[i], p_opens[i + 1], "P_open must be monotonically increasing with dose")

        for po in p_opens:
            self.assertGreaterEqual(po, 0.0)
            self.assertLessEqual(po, 1.0)

    def test_11_dose_response_surface_generation(self):
        """compute_dose_response_surface must generate a grid of valid P_open points."""
        surface = compute_dose_response_surface(
            variant_id="F508del",
            compound_x="VX-445",
            compound_y="VX-770",
            conc_range_x=[0.1, 0.5, 1.0],
            conc_range_y=[0.05, 0.2, 0.5],
        )
        self.assertEqual(len(surface), 9)
        for pt in surface:
            self.assertIn("conc_x", pt)
            self.assertIn("conc_y", pt)
            self.assertIn("p_open", pt)
            self.assertGreaterEqual(pt["p_open"], 0.0)

    def test_12_disclaimers_present_in_results(self):
        """All output records must include NLM and Medical disclaimers."""
        cycle = AllostericModulatorCycle(variant_id="F508del")
        res = cycle.run_cycle(["VX-770"])
        self.assertEqual(res.nl_registry_disclaimer, NLM_REGISTRY_DISCLAIMER)
        self.assertEqual(res.medical_disclaimer, MEDICAL_DISCLAIMER)

    def test_13_cli_execution_and_json_mode(self):
        """Script must execute via CLI with both human-readable and JSON output formats."""
        script_path = REPO_ROOT / "tools" / "discovery" / "allosteric_potentiator_model.py"
        cmd = [sys.executable, str(script_path), "--variant", "F508del", "--regimen", "Trikafta", "--json"]
        out = subprocess.check_output(cmd, universal_newlines=True)
        data = json.loads(out)

        self.assertEqual(data["variant_id"], "F508del")
        self.assertEqual(data["modulators_applied"], ["VX-445", "VX-661", "VX-770"])
        self.assertIn("thermodynamic_synergy_s", data)
        self.assertIn("predicted_delta_swcl_mmol_L", data)


if __name__ == "__main__":
    unittest.main()
