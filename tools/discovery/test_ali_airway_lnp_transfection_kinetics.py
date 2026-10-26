#!/usr/bin/env python3
"""Unit test suite for ALI airway epithelial aerosolized LNP transfection kinetics (T-0446)."""
from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ali_airway_lnp_transfection_kinetics import (
    ALIAirwayParameters,
    CLINICAL_TRIAL_4D710,
    CLINICAL_TRIAL_VX522,
    NLM_ATTRIBUTION,
    PDB_ACTIVE_CFTR,
    PDB_DEPHOSPHO_CFTR,
    PRIMARY_DOI_AIRWAY_PHYSIOLOGY,
    PRIMARY_DOI_INHALED_LNP,
    PRIMARY_DOI_MUCUS_BARRIER,
    PRIMARY_DOI_MUCUS_DIFFUSION,
    calculate_effective_mucus_diffusivity,
    calculate_electrostatic_penalty,
    calculate_mucus_steric_diffusivity,
    calculate_renkin_hindrance_partitioning,
    calculate_stokes_einstein_diffusion,
    evaluate_tight_junction_permeation,
    simulate_cftr_expression_and_isc_rescue,
    simulate_full_ali_airway_transfection,
    simulate_mucociliary_filtration_competition,
    simulate_polarized_cellular_internalization,
)


class TestALIAirwayLNPTransfectionKinetics(unittest.TestCase):
    """Test suite verifying biophysical accuracy, grounding, and clinical safety boundaries."""

    def test_stokes_einstein_and_ogston_sieving(self) -> None:
        """Verify unhindered Stokes-Einstein diffusion and steric mesh filtration."""
        # 65 nm LNP in water at 37 C
        d0 = calculate_stokes_einstein_diffusion(65.0)
        # Should be approximately 1.0e-11 m^2/s
        self.assertGreater(d0, 5.0e-12)
        self.assertLess(d0, 2.0e-11)

        # Negative diameter raises ValueError
        with self.assertRaises(ValueError):
            calculate_stokes_einstein_diffusion(-10.0)

        # Steric hindrance through 120 nm pore mesh
        s_small = calculate_mucus_steric_diffusivity(30.0, mesh_size_nm=120.0)
        s_med = calculate_mucus_steric_diffusivity(65.0, mesh_size_nm=120.0)
        s_large = calculate_mucus_steric_diffusivity(250.0, mesh_size_nm=120.0)

        self.assertGreater(s_small, s_med)
        self.assertGreater(s_med, 0.0)
        self.assertEqual(s_large, 0.0)

    def test_electrostatic_penalty(self) -> None:
        """Verify neutral / mildly anionic LNPs evade mucin charge binding while cationic particles are trapped."""
        e_neutral = calculate_electrostatic_penalty(-1.5)
        e_cationic = calculate_electrostatic_penalty(+12.0)

        self.assertGreater(e_neutral, 0.90)
        self.assertLess(e_cationic, 0.10)

    def test_mucociliary_shear_filtration_sweep(self) -> None:
        """Verify velocity sweep (5-20 mm/min) maintains productive retention."""
        for v_mcc in [5.0, 10.0, 15.0, 20.0]:
            res = simulate_mucociliary_filtration_competition(
                diameter_nm=65.0,
                zeta_potential_mv=-1.5,
                mucus_thickness_um=25.0,
                mcc_velocity_mm_per_min=v_mcc,
                airway_path_length_mm=15.0,
            )
            retention = res["productive_epithelial_retention_fraction"]
            self.assertGreater(retention, 0.0)
            self.assertLessEqual(retention, 1.0)
            # Faster clearance slightly reduces mucosal retention
            if v_mcc == 5.0:
                self.assertGreater(retention, 0.15)

    def test_endosomal_escape_threshold_guard(self) -> None:
        """Verify model enforces > 5.0% endosomal escape efficiency."""
        # 8.0% escape passes
        res = simulate_polarized_cellular_internalization(
            retention_fraction=0.30,
            dose_mrna_copies_per_cm2=1.0e8,
            endosomal_escape_efficiency=0.08,
        )
        self.assertGreater(res["cytosolic_functional_mrna_copies"], 0.0)
        self.assertEqual(res["endosomal_escape_efficiency"], 0.08)

        # <= 5.0% fails-closed with ValueError
        with self.assertRaises(ValueError):
            simulate_polarized_cellular_internalization(
                retention_fraction=0.30,
                dose_mrna_copies_per_cm2=1.0e8,
                endosomal_escape_efficiency=0.05,
            )

    def test_tight_junction_dynamics_and_teer_safety(self) -> None:
        """Verify transient ZO-1/claudin-1 opening enables basal access while preserving TER >= 500 ohm*cm^2."""
        tj = evaluate_tight_junction_permeation(
            lnp_diameter_nm=65.0,
            dilated_pore_radius_nm=35.0,
            resting_pore_radius_nm=0.5,
            baseline_teer_ohm_cm2=1200.0,
            min_safe_teer_ohm_cm2=500.0,
        )
        self.assertTrue(tj["barrier_integrity_preserved"])
        self.assertGreaterEqual(tj["effective_teer_ohm_cm2"], 500.0)
        self.assertGreater(tj["renkin_partitioning_factor"], 0.0)
        self.assertGreater(tj["paracellular_basal_permeability_ratio"], 0.0)

    def test_cftr_ribosomal_translation_and_isc_rescue(self) -> None:
        """Verify CFTR translation yields therapeutic functional rescue >= 10% wild-type."""
        rescue = simulate_cftr_expression_and_isc_rescue(
            cytosolic_mrna_copies_per_cm2=5.0e6,
            translation_rate_codons_s=4.2,
            cftr_codons=1480,
            folding_efficiency=0.32,
            baseline_isc_ua_cm2=1.0,
            wt_isc_ua_cm2=30.0,
        )
        self.assertTrue(rescue["cftr_rescue_adequate_ge_10pct"])
        self.assertGreaterEqual(rescue["percent_wildtype_rescue"], 10.0)
        self.assertGreater(rescue["rescued_delta_isc_ua_cm2"], 2.9)
        self.assertAlmostEqual(rescue["translation_duration_seconds"], 1480.0 / 4.2, places=1)

    def test_end_to_end_ali_transfection_simulation(self) -> None:
        """Verify full end-to-end transfection simulation satisfies all acceptance invariants."""
        params = ALIAirwayParameters(
            lnp_diameter_nm=65.0,
            zeta_potential_mv=-1.5,
            mcc_velocity_mm_per_min=10.0,
            endosomal_escape_efficiency=0.08,
        )
        sim = simulate_full_ali_airway_transfection(params, dose_mrna_copies_per_cm2=1.0e9)

        self.assertEqual(sim["status"], "PASS")
        self.assertIn(CLINICAL_TRIAL_VX522, sim["clinical_trial_grounding"])
        self.assertIn(CLINICAL_TRIAL_4D710, sim["clinical_trial_grounding"])
        self.assertIn(PDB_DEPHOSPHO_CFTR, sim["pdb_structures"])
        self.assertIn(PDB_ACTIVE_CFTR, sim["pdb_structures"])
        self.assertIn(PRIMARY_DOI_MUCUS_DIFFUSION, sim["literature_citations"])
        self.assertIn(PRIMARY_DOI_AIRWAY_PHYSIOLOGY, sim["literature_citations"])
        self.assertIn(PRIMARY_DOI_MUCUS_BARRIER, sim["literature_citations"])
        self.assertIn(PRIMARY_DOI_INHALED_LNP, sim["literature_citations"])
        self.assertIn("ClinicalTrials.gov", sim["nlm_attribution"])

        # Check sub-module assertions
        self.assertTrue(sim["tight_junction_paracellular_access"]["barrier_integrity_preserved"])
        self.assertGreaterEqual(
            sim["tight_junction_paracellular_access"]["effective_teer_ohm_cm2"], 500.0
        )
        self.assertGreater(
            sim["polarized_internalization"]["endosomal_escape_efficiency"], 0.05
        )
        self.assertTrue(
            sim["functional_cftr_rescue"]["cftr_rescue_adequate_ge_10pct"]
        )
        self.assertGreaterEqual(
            sim["functional_cftr_rescue"]["percent_wildtype_rescue"], 10.0
        )


if __name__ == "__main__":
    unittest.main()
