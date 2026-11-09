#!/usr/bin/env python3
"""Unit tests for Round 153: Non-Equilibrium Brinkman-Stokes Nanoparticle Mobility & Mucus Steric/Electrostatic Filtration Invariant Bridge (T-0810)."""
from __future__ import annotations

import unittest

from tools.discovery.round153_mucus_nanoparticle_filtration_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    compute_brinkman_drag_factor,
    compute_electrostatic_escape_factor,
    compute_ogston_sieving_factor,
    compute_stokes_einstein_diffusivity,
    evaluate_falsification_negative_control,
    evaluate_round153_mucus_nanoparticle_filtration,
)


class Round153NanoparticleFiltrationTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round153_mucus_nanoparticle_filtration(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.free_diffusivity_um2_s, 0.0)
            self.assertGreater(res.macroscopic_diffusivity_um2_s, 0.0)
            self.assertGreater(res.ogston_steric_sieving_factor, 0.0)
            self.assertLessEqual(res.ogston_steric_sieving_factor, 1.0)
            self.assertGreater(res.brinkman_porous_drag_factor, 0.0)
            self.assertLessEqual(res.brinkman_porous_drag_factor, 1.0)
            self.assertGreater(res.electrostatic_escape_factor, 0.0)
            self.assertLessEqual(res.electrostatic_escape_factor, 1.0)
            self.assertGreaterEqual(res.first_passage_penetration_probability, 0.0)
            self.assertLessEqual(res.first_passage_penetration_probability, 1.0)

    def test_healthy_airway_stealth_lnp(self) -> None:
        healthy = evaluate_round153_mucus_nanoparticle_filtration("healthy_airway_stealth_lnp")
        self.assertTrue(healthy.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(healthy.macroscopic_diffusivity_um2_s, 3.5)
        self.assertGreaterEqual(healthy.first_passage_penetration_probability, 0.85)
        self.assertAlmostEqual(healthy.fractional_wild_type_penetration, 1.0, delta=0.05)
        self.assertLess(healthy.diffusion_time_to_pcl_s, 60.0)

    def test_cf_entrapment_modes(self) -> None:
        cf_stealth = evaluate_round153_mucus_nanoparticle_filtration("cf_hyperconcentrated_mucus_stealth_lnp")
        cf_cationic = evaluate_round153_mucus_nanoparticle_filtration("cf_cationic_lipoplex_entrapment")

        # Stealth in CF has reduced mobility due to pore constriction
        self.assertLess(cf_stealth.macroscopic_diffusivity_um2_s, 1.5)
        self.assertGreater(cf_stealth.macroscopic_diffusivity_um2_s, 0.5)

        # Cationic lipoplex is completely trapped
        self.assertFalse(cf_cationic.crosses_kahn_kalai_threshold)
        self.assertLess(cf_cationic.macroscopic_diffusivity_um2_s, 0.01)
        self.assertLess(cf_cationic.first_passage_penetration_probability, 0.01)
        self.assertGreater(cf_cationic.diffusion_time_to_pcl_s, 1e5)

    def test_kahn_kalai_filtration_boundary(self) -> None:
        kk = evaluate_round153_mucus_nanoparticle_filtration("kahn_kalai_filtration_boundary")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.first_passage_penetration_probability, 0.118)
        self.assertGreaterEqual(kk.macroscopic_diffusivity_um2_s, 0.50)

    def test_mucolytic_dornase_alfa_rescue(self) -> None:
        rescued = evaluate_round153_mucus_nanoparticle_filtration("mucolytic_dornase_alfa_rescue")
        self.assertTrue(rescued.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(rescued.macroscopic_diffusivity_um2_s, 3.0)
        self.assertGreaterEqual(rescued.first_passage_penetration_probability, 0.80)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_particle_radius=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")

        poisoned_r = evaluate_falsification_negative_control(poison_negative_particle_radius=True)
        self.assertFalse(poisoned_r["passed"])

        poisoned_xi = evaluate_falsification_negative_control(poison_zero_pore_spacing=True)
        self.assertFalse(poisoned_xi["passed"])

    def test_scaling_and_monotonicity_laws(self) -> None:
        # Stokes-Einstein size scaling
        d0_30nm = compute_stokes_einstein_diffusivity(15e-9)
        d0_65nm = compute_stokes_einstein_diffusivity(32.5e-9)
        d0_135nm = compute_stokes_einstein_diffusivity(67.5e-9)
        self.assertGreater(d0_30nm, d0_65nm)
        self.assertGreater(d0_65nm, d0_135nm)

        # Ogston pore sieving
        ogston_wide = compute_ogston_sieving_factor(32.5, 200.0)
        ogston_tight = compute_ogston_sieving_factor(32.5, 70.0)
        self.assertGreater(ogston_wide, ogston_tight)

        # Brinkman drag
        kappa_wide = (200e-9 ** 2) / 16.0
        kappa_tight = (70e-9 ** 2) / 16.0
        brinkman_wide = compute_brinkman_drag_factor(32.5e-9, kappa_wide)
        brinkman_tight = compute_brinkman_drag_factor(32.5e-9, kappa_tight)
        self.assertGreater(brinkman_wide, brinkman_tight)

        # Electrostatic escape
        f_stealth = compute_electrostatic_escape_factor(-0.8)
        f_cationic = compute_electrostatic_escape_factor(+14.5)
        self.assertGreater(f_stealth, f_cationic)
        self.assertLess(f_cationic, 0.01)

    def test_two_column_epistemic_hygiene_and_five_layer_ladder(self) -> None:
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("falsification_condition", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        for layer in (
            "cool_human_fact",
            "everyday_analogy",
            "engineering_lens",
            "formal_theorem",
            "falsifiable_code_and_toy",
        ):
            self.assertIn(layer, FIVE_LAYER_LADDER)
            self.assertGreaterEqual(len(FIVE_LAYER_LADDER[layer]), 25)
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 4)


if __name__ == "__main__":
    unittest.main()
