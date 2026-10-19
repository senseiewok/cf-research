#!/usr/bin/env python3
"""Deterministic verification test suite for aav_tropism_model.py (T-0208).

Dual test runner compatibility:
- pytest tools/discovery/test_aav_tropism_model.py -v
- python -m unittest discover -s tools/discovery -p "test_*.py"
"""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import aav_tropism_model as atm


class TestAerosolDeposition(unittest.TestCase):
    """Verifies aerodynamic droplet regional lung deposition mechanics."""

    def test_deposition_mass_conservation(self):
        """All regional deposition fractions must sum exactly to 1.0."""
        for mmad in [1.5, 2.5, 3.5, 4.5, 6.0, 8.0]:
            dep = atm.calculate_aerosol_deposition(mmad_um=mmad, gsd=1.8, inhalation_flow_l_min=15.0)
            total = (
                dep.extrathoracic_fraction
                + dep.tracheobronchial_fraction
                + dep.alveolar_fraction
                + dep.exhaled_fraction
            )
            self.assertAlmostEqual(total, 1.0, places=3)

    def test_extrathoracic_impaction_monotonicity(self):
        """Larger droplets must exhibit higher extrathoracic (oropharyngeal) impaction."""
        dep_small = atm.calculate_aerosol_deposition(mmad_um=2.0)
        dep_large = atm.calculate_aerosol_deposition(mmad_um=6.0)
        self.assertGreater(dep_large.extrathoracic_fraction, dep_small.extrathoracic_fraction)

    def test_tracheobronchial_window(self):
        """Conducting airway deposition must be substantial in the 3.0-4.0 um clinical nebulizer window."""
        dep_opt = atm.calculate_aerosol_deposition(mmad_um=3.5)
        self.assertGreaterEqual(dep_opt.tracheobronchial_fraction, 0.15)

    def test_alveolar_fine_droplet_dominance(self):
        """Alveolar deposition should be higher for fine droplets (2.0 um) than coarse droplets (7.0 um)."""
        dep_fine = atm.calculate_aerosol_deposition(mmad_um=2.0)
        dep_coarse = atm.calculate_aerosol_deposition(mmad_um=7.0)
        self.assertGreater(dep_fine.alveolar_fraction, dep_coarse.alveolar_fraction)

    def test_dose_distribution_partition(self):
        """Nominal dose must partition correctly through nebulizer efficiency."""
        dist = atm.calculate_aerosol_dose_distribution(
            nominal_dose_vg=1.0e13,
            nebulizer_efficiency=0.60,
            mmad_um=3.5,
        )
        self.assertAlmostEqual(dist.inhaled_dose_vg, 6.0e12, delta=1.0e9)
        total_delivered = dist.extrathoracic_vg + dist.tracheobronchial_vg + dist.alveolar_vg + dist.exhaled_vg
        self.assertAlmostEqual(total_delivered, dist.inhaled_dose_vg, delta=1.0e9)

    def test_invalid_aerosol_parameters(self):
        """Must reject non-physical inputs."""
        with self.assertRaises(ValueError):
            atm.calculate_aerosol_deposition(mmad_um=-2.0)
        with self.assertRaises(ValueError):
            atm.calculate_aerosol_deposition(gsd=0.8)  # GSD must be >= 1.0
        with self.assertRaises(ValueError):
            atm.calculate_aerosol_deposition(inhalation_flow_l_min=-10.0)
        with self.assertRaises(ValueError):
            atm.calculate_aerosol_dose_distribution(nominal_dose_vg=-1.0)
        with self.assertRaises(ValueError):
            atm.calculate_aerosol_dose_distribution(nominal_dose_vg=1.0e12, nebulizer_efficiency=1.5)


class TestMucusHindranceAndDiffusion(unittest.TestCase):
    """Verifies mucin mesh pore scaling, Stokes-Einstein diffusion, and MCC traversal."""

    def test_mesh_pore_deconcentration_scaling(self):
        """Mucus pore size must contract monotonically with increasing solids (Hill et al.)."""
        prof_healthy = atm.calculate_mucus_diffusion_profile(solids_pct=2.0)
        prof_cf = atm.calculate_mucus_diffusion_profile(solids_pct=7.5)
        prof_severe = atm.calculate_mucus_diffusion_profile(solids_pct=10.0)

        self.assertGreater(prof_healthy.mesh_pore_diameter_nm, prof_cf.mesh_pore_diameter_nm)
        self.assertGreater(prof_cf.mesh_pore_diameter_nm, prof_severe.mesh_pore_diameter_nm)
        self.assertAlmostEqual(prof_healthy.mesh_pore_diameter_nm, 180.0, delta=1.0)
        # CF 7.5% solids should compress pores to under 70 nm
        self.assertLess(prof_cf.mesh_pore_diameter_nm, 70.0)

    def test_capsid_mucin_adhesion_impact(self):
        """Directed evolution capsid (low adhesion) must diffuse faster than wild-type AAV2 (high adhesion)."""
        prof_4d = atm.calculate_mucus_diffusion_profile(
            solids_pct=7.5,
            capsid_adhesion_coeff=atm.CAPSID_CATALOG["4D-A101"].mucin_adhesion_coefficient,
        )
        prof_aav2 = atm.calculate_mucus_diffusion_profile(
            solids_pct=7.5,
            capsid_adhesion_coeff=atm.CAPSID_CATALOG["AAV2"].mucin_adhesion_coefficient,
        )
        # 4D-A101 effective diffusivity should exceed AAV2 by at least 4-fold in CF mucus
        self.assertGreater(
            prof_4d.effective_diffusivity_um2_s,
            prof_aav2.effective_diffusivity_um2_s * 4.0,
        )
        # Mucus penetration fraction must also be substantially higher
        self.assertGreater(
            prof_4d.mucus_penetration_fraction,
            prof_aav2.mucus_penetration_fraction,
        )

    def test_mucus_thickness_timescale_scaling(self):
        """Diffusion timescale must scale quadratically with mucus layer thickness (tau propto L^2)."""
        prof_thin = atm.calculate_mucus_diffusion_profile(mucus_thickness_um=15.0)
        prof_thick = atm.calculate_mucus_diffusion_profile(mucus_thickness_um=30.0)
        ratio = prof_thick.diffusion_timescale_s / prof_thin.diffusion_timescale_s
        self.assertAlmostEqual(ratio, 4.0, places=1)

    def test_invalid_mucus_inputs(self):
        with self.assertRaises(ValueError):
            atm.calculate_mucus_diffusion_profile(solids_pct=0.0)
        with self.assertRaises(ValueError):
            atm.calculate_mucus_diffusion_profile(solids_pct=25.0)
        with self.assertRaises(ValueError):
            atm.calculate_mucus_diffusion_profile(capsid_adhesion_coeff=-0.5)
        with self.assertRaises(ValueError):
            atm.calculate_mucus_diffusion_profile(mucus_thickness_um=-5.0)


class TestNeutralizingAntibodyEvasion(unittest.TestCase):
    """Verifies mucosal antibody transudation and log-logistic neutralization kinetics."""

    def test_zero_titer_full_activity(self):
        """Zero serum neutralizing antibody titer must yield 100% active capsids."""
        nab = atm.calculate_neutralizing_antibody_profile(
            serum_nab_titer=0.0,
            capsid_ic50_serum_titer=50.0,
        )
        self.assertEqual(nab.active_unneutralized_fraction, 1.0)
        self.assertEqual(nab.neutralization_fraction, 0.0)

    def test_titration_monotonic_neutralization(self):
        """Increasing antibody titer must monotonically decrease active unneutralized capsid fraction."""
        titers = [1.0, 10.0, 50.0, 200.0, 1000.0]
        active_fractions = [
            atm.calculate_neutralizing_antibody_profile(
                serum_nab_titer=t,
                capsid_ic50_serum_titer=50.0,
            ).active_unneutralized_fraction
            for t in titers
        ]
        for i in range(len(active_fractions) - 1):
            self.assertGreater(active_fractions[i], active_fractions[i + 1])

    def test_4da101_nab_evasion_advantage(self):
        """4D-A101 must maintain higher active fraction than AAV2 in presence of moderate antibodies (titer 1:20)."""
        nab_4d = atm.calculate_neutralizing_antibody_profile(
            serum_nab_titer=20.0,
            capsid_ic50_serum_titer=atm.CAPSID_CATALOG["4D-A101"].ic50_serum_nab_titer,
            hill_coefficient=atm.CAPSID_CATALOG["4D-A101"].hill_coefficient,
        )
        nab_aav2 = atm.calculate_neutralizing_antibody_profile(
            serum_nab_titer=20.0,
            capsid_ic50_serum_titer=atm.CAPSID_CATALOG["AAV2"].ic50_serum_nab_titer,
            hill_coefficient=atm.CAPSID_CATALOG["AAV2"].hill_coefficient,
        )
        self.assertGreater(nab_4d.active_unneutralized_fraction, nab_aav2.active_unneutralized_fraction)
        self.assertGreater(nab_4d.active_unneutralized_fraction, 0.80)

    def test_invalid_nab_inputs(self):
        with self.assertRaises(ValueError):
            atm.calculate_neutralizing_antibody_profile(serum_nab_titer=-5.0, capsid_ic50_serum_titer=50.0)
        with self.assertRaises(ValueError):
            atm.calculate_neutralizing_antibody_profile(serum_nab_titer=10.0, capsid_ic50_serum_titer=-1.0)
        with self.assertRaises(ValueError):
            atm.calculate_neutralizing_antibody_profile(serum_nab_titer=10.0, capsid_ic50_serum_titer=50.0, hill_coefficient=-0.5)
        with self.assertRaises(ValueError):
            atm.calculate_neutralizing_antibody_profile(serum_nab_titer=10.0, capsid_ic50_serum_titer=50.0, mucosal_transudation_ratio=1.5)


class TestTransductionAndCellularPoisson(unittest.TestCase):
    """Verifies end-to-end airway transduction and Poisson cellular transgene coverage."""

    def test_4da101_clinical_regime_outcome(self):
        """Evaluates 4D-A101 under typical clinical trial dose parameters (NCT05248230)."""
        capsid = atm.CAPSID_CATALOG["4D-A101"]
        outcome = atm.calculate_transduction_outcome(
            capsid=capsid,
            nominal_dose_vg=1.0e13,
            serum_nab_titer=20.0,
            mucus_solids_pct=7.5,
        )
        self.assertEqual(outcome.capsid_name, "4D-A101")
        self.assertGreater(outcome.effective_moi_vg_per_cell, 0.05)
        self.assertGreater(outcome.transduced_cells_fraction_poisson, 0.05)
        self.assertGreater(outcome.estimated_transduced_cell_count, 1.0e7)
        self.assertIn("ClinicalTrials.gov does not evaluate safety or efficacy", outcome.nlm_disclaimer)

    def test_dose_dependent_poisson_saturation(self):
        """Increasing nominal dose must monotonically increase Poisson cellular coverage toward 1.0."""
        capsid = atm.CAPSID_CATALOG["4D-A101"]
        doses = [1.0e10, 1.0e11, 1.0e12, 1.0e13]
        coverages = [
            atm.calculate_transduction_outcome(capsid=capsid, nominal_dose_vg=d).transduced_cells_fraction_poisson
            for d in doses
        ]
        for i in range(len(coverages) - 1):
            self.assertGreater(coverages[i + 1], coverages[i])
        self.assertGreaterEqual(coverages[-1], 0.80)

    def test_4da101_superiority_over_wt_aav2(self):
        """In CF mucus with antibodies, 4D-A101 must achieve substantially higher transduction than AAV2."""
        outcome_4d = atm.calculate_transduction_outcome(
            capsid=atm.CAPSID_CATALOG["4D-A101"],
            nominal_dose_vg=1.0e13,
            serum_nab_titer=50.0,
            mucus_solids_pct=7.5,
        )
        outcome_aav2 = atm.calculate_transduction_outcome(
            capsid=atm.CAPSID_CATALOG["AAV2"],
            nominal_dose_vg=1.0e13,
            serum_nab_titer=50.0,
            mucus_solids_pct=7.5,
        )
        self.assertGreater(
            outcome_4d.transduced_cells_fraction_poisson,
            outcome_aav2.transduced_cells_fraction_poisson * 5.0,
        )


class TestCapsidCatalogAndCLI(unittest.TestCase):
    """Verifies catalog integrity, comparison matrix, and CLI execution."""

    def test_catalog_completeness(self):
        """All catalog entries must have valid biological parameters and descriptions."""
        self.assertIn("4D-A101", atm.CAPSID_CATALOG)
        self.assertIn("AAV2", atm.CAPSID_CATALOG)
        self.assertIn("AAV5", atm.CAPSID_CATALOG)
        self.assertIn("AAV6", atm.CAPSID_CATALOG)
        self.assertIn("AAV9", atm.CAPSID_CATALOG)
        self.assertIn("AAV2.5T", atm.CAPSID_CATALOG)

        for name, prop in atm.CAPSID_CATALOG.items():
            self.assertEqual(prop.name, name)
            self.assertGreater(prop.apical_receptor_tropism_score, 0.0)
            self.assertGreater(prop.ic50_serum_nab_titer, 0.0)
            self.assertGreater(prop.intracellular_trafficking_efficiency, 0.0)

    def test_comparative_matrix_generation(self):
        """compare_capsids_across_titer_spectrum must evaluate all requested capsids."""
        matrix = atm.compare_capsids_across_titer_spectrum(
            nominal_dose_vg=5.0e12,
            serum_titers=[1.0, 20.0, 100.0],
            capsid_names=["4D-A101", "AAV2"],
        )
        self.assertIn("capsids", matrix)
        self.assertIn("4D-A101", matrix["capsids"])
        self.assertIn("AAV2", matrix["capsids"])
        self.assertEqual(len(matrix["capsids"]["4D-A101"]["titer_evaluations"]), 3)

    def test_cli_execution_json(self):
        """Executes CLI with --json and asserts valid structured JSON stdout."""
        script_path = HERE / "aav_tropism_model.py"
        res = subprocess.run(
            [sys.executable, str(script_path), "--capsid", "4D-A101", "--json"],
            capture_output=True,
            text=True,
            check=True,
        )
        data = json.loads(res.stdout)
        self.assertEqual(data["capsid_name"], "4D-A101")
        self.assertIn("effective_moi_vg_per_cell", data)
        self.assertIn("nlm_disclaimer", data)

    def test_cli_execution_compare_all(self):
        """Executes CLI with --compare-all --json."""
        script_path = HERE / "aav_tropism_model.py"
        res = subprocess.run(
            [sys.executable, str(script_path), "--compare-all", "--json"],
            capture_output=True,
            text=True,
            check=True,
        )
        data = json.loads(res.stdout)
        self.assertIn("capsids", data)
        self.assertIn("4D-A101", data["capsids"])


class TestAntiHallucinationAndProvenance(unittest.TestCase):
    """Verifies strict grounding against clinical registries, primary DOIs, and disclaimer rules."""

    def test_nct05248230_in_snapshot(self):
        """Clinical trial NCT05248230 must exist in trial_atlas/snapshot.json."""
        snapshot_file = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
        self.assertTrue(snapshot_file.exists(), f"Missing {snapshot_file}")
        snap = json.loads(snapshot_file.read_text(encoding="utf-8"))
        snapshot_ncts = {s["nct_id"] for s in snap.get("studies", []) if "nct_id" in s}
        self.assertIn("NCT05248230", snapshot_ncts)

    def test_nlm_disclaimer_verbatim(self):
        """NLM disclaimer must match exact legal formulation."""
        self.assertIn(
            "ClinicalTrials.gov does not evaluate safety or efficacy",
            atm.NLM_REGISTRY_DISCLAIMER,
        )
        self.assertIn(
            "listing on the registry does not denote government endorsement",
            atm.NLM_REGISTRY_DISCLAIMER,
        )

    def test_no_prohibited_overclaims(self):
        """Module source code must contain zero prohibited clinical overclaim phrases."""
        script_path = HERE / "aav_tropism_model.py"
        content = script_path.read_text(encoding="utf-8").lower()
        prohibited = [
            "proven to cure",
            "guarantees clinical benefit",
            "recommended treatment regimen",
            "safe and effective for patients",
            "cures cystic fibrosis",
        ]
        for phrase in prohibited:
            self.assertNotIn(phrase, content, f"Detected prohibited overclaim: '{phrase}'")


if __name__ == "__main__":
    unittest.main()
