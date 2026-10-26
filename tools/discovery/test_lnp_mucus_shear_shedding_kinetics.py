#!/usr/bin/env python3
"""Unit tests for LNP Mucus Fluidic Shear & PEG Shedding Kinetics (T-0439)."""
import math
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lnp_mucus_shear_shedding_kinetics import (
    CLINICAL_TRIAL_4D710,
    CLINICAL_TRIAL_VX522,
    NLM_ATTRIBUTION,
    PDB_ACTIVE_CFTR,
    PDB_DEPHOSPHO_CFTR,
    NebulizerShearParameters,
    PEGSheddingKinetics,
    calculate_dlvo_electrostatic_factor,
    calculate_membrane_tension,
    calculate_nebulizer_shear_stress,
    calculate_ogston_renkin_diffusivity,
    evaluate_lnp_nebulization_survival,
    simulate_mucus_transit_and_basal_delivery,
)


class TestLNPMucusShearSheddingKinetics(unittest.TestCase):
    """Test suite verifying nebulizer shear stability, steric sieving, DLVO, and shedding kinetics."""

    def test_nebulizer_shear_stress_and_tension(self) -> None:
        """Verify shear stress and membrane tension scale linearly with shear rate and radius."""
        shear_rate = 1.0e5  # 10^5 s^-1
        tau = calculate_nebulizer_shear_stress(shear_rate, viscosity_pa_s=0.001)
        self.assertAlmostEqual(tau, 100.0, places=2)  # 100 Pa

        radius_nm = 32.5  # 65 nm diameter
        tension = calculate_membrane_tension(radius_nm, tau)
        # 0.5 * 32.5e-9 * 100 * 1000 = 1.625 mN/m
        self.assertAlmostEqual(tension, 1.625, places=3)

    def test_nebulization_survival_regimes(self) -> None:
        """Verify LNP remains intact under standard nebulization but lyses under extreme shear."""
        # Standard vibrating mesh shear: 10^5 s^-1 -> intact
        res_standard = evaluate_lnp_nebulization_survival(32.5, 1.0e5)
        self.assertEqual(res_standard["verdict"], "INTACT")
        self.assertGreaterEqual(res_standard["survival_fraction"], 0.95)

        # Extreme non-physiological shear: 5.0 * 10^6 s^-1 -> membrane tension > critical lysis (20 mN/m) -> lyses
        res_extreme = evaluate_lnp_nebulization_survival(32.5, 5.0e6)
        self.assertEqual(res_extreme["verdict"], "LYSED")
        self.assertLess(res_extreme["survival_fraction"], 0.50)

    def test_ogston_renkin_steric_diffusivity_bounds(self) -> None:
        """Verify steric diffusivity ratio is monotonic and drops to zero at cutoff."""
        mesh_nm = 120.0
        # Small particle (20 nm) diffuses well
        d_small = calculate_ogston_renkin_diffusivity(20.0, mesh_size_nm=mesh_nm)
        # Moderate particle (65 nm) has reduced diffusivity
        d_med = calculate_ogston_renkin_diffusivity(65.0, mesh_size_nm=mesh_nm)
        # Large particle (250 nm > 2 * mesh_nm) completely excluded
        d_large = calculate_ogston_renkin_diffusivity(250.0, mesh_size_nm=mesh_nm)

        self.assertGreater(d_small, d_med)
        self.assertGreater(d_med, 0.0)
        self.assertEqual(d_large, 0.0)

    def test_dlvo_electrostatic_neutral_vs_cationic(self) -> None:
        """Verify near-neutral PEGylated LNPs pass freely while cationic particles are entrapped."""
        # Near-neutral (-1.0 mV)
        e_neutral = calculate_dlvo_electrostatic_factor(-1.0)
        self.assertGreater(e_neutral, 0.85)

        # Strongly cationic (+15.0 mV) binds sialylated mucin
        e_cationic = calculate_dlvo_electrostatic_factor(+15.0)
        self.assertLess(e_cationic, 0.05)

    def test_peg_shedding_kinetics(self) -> None:
        """Verify acyl-chain shedding kinetics enable endocytosis competence over time."""
        # DMG-PEG (k = 0.30/h) sheds from 2.0% down to <= 0.5% in ~5 hours
        dmg = PEGSheddingKinetics(anchor_type="DMG_PEG2000", k_shed_per_hour=0.30)
        self.assertFalse(dmg.is_endocytosis_competent(0.0))
        # At t = 5 hours, 2.0 * exp(-1.5) = 0.446% <= 0.5%
        self.assertTrue(dmg.is_endocytosis_competent(5.0))

    def test_coupled_mucus_transit_and_basal_delivery(self) -> None:
        """Verify end-to-end delivery of 65nm near-neutral LNP to the basal niche."""
        sim = simulate_mucus_transit_and_basal_delivery(
            diameter_lnp_nm=65.0,
            zeta_potential_mv=-1.2,
            anchor_type="DMG_PEG2000",
            k_shed_per_hour=0.30,
        )
        self.assertTrue(sim["curative_delivery_feasible"])
        self.assertGreater(sim["effective_diffusivity_ratio"], 0.10)
        self.assertGreater(sim["delivered_fraction_at_basal"], 0.10)


if __name__ == "__main__":
    unittest.main()
