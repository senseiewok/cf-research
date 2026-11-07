"""Unit tests for Round 108: Non-Equilibrium Kramers Escape & Stokes-Brinkman Viscoelastic Mucus Mesh Gel Filtration Dynamics."""

import unittest
from tools.discovery.round108_kramers_stokes_brinkman_mucus_filtration import (
    CITATIONS,
    HYPOTHESIS,
    MucusFiltrationProfile,
    PARAMETER_PROVENANCE,
    PROFILES,
    ROUND_ID,
    TASK_ID,
    TITLE,
    run_falsification,
    run_falsification_checks,
    simulate_mucus_filtration,
)


class TestRound108KramersStokesBrinkmanMucusFiltration(unittest.TestCase):
    def test_01_metadata_and_ah10_contract(self) -> None:
        self.assertEqual(ROUND_ID, "round108")
        self.assertEqual(TASK_ID, "T-0737")
        self.assertIn("Kramers Escape", TITLE)
        self.assertIn("Stokes-Brinkman", TITLE)
        self.assertIn("dornase alfa", HYPOTHESIS)
        self.assertEqual(len(CITATIONS), 5)
        self.assertTrue(all(c.get("scite_verified") is True for c in CITATIONS))

        ah10 = run_falsification()
        self.assertEqual(ah10["rule"], "AH-10")
        self.assertEqual(ah10["round_id"], "round108")
        self.assertEqual(ah10["task_id"], "T-0737")
        self.assertTrue(ah10["passed"])

    def test_02_parameter_provenance_completeness(self) -> None:
        required_keys = (
            "particle_radius_nm",
            "mesh_pore_size_nm",
            "brinkman_permeability_nm2",
            "adhesive_barrier_kbt",
            "mucus_thickness_um",
        )
        for key in required_keys:
            self.assertIn(key, PARAMETER_PROVENANCE)
            meta = PARAMETER_PROVENANCE[key]
            self.assertIn("unit", meta)
            self.assertIn("bounds", meta)
            self.assertIn("source", meta)

    def test_03_falsification_checks_all_pass(self) -> None:
        res = run_falsification_checks()
        self.assertEqual(res["round"], 108)
        self.assertTrue(res["passed"])
        self.assertEqual(len(res["results"]), 5)
        for name, ok in res["checks"].items():
            self.assertTrue(ok, f"Failed check: {name}")

    def test_04_stokes_brinkman_and_ogston_physical_bounds(self) -> None:
        for profile in PROFILES:
            sim = simulate_mucus_filtration(profile)
            self.assertTrue(sim.is_stokes_brinkman_bounded)
            self.assertTrue(sim.is_kramers_barrier_passive)
            self.assertGreater(sim.ogston_steric_factor, 0.0)
            self.assertLessEqual(sim.ogston_steric_factor, 1.0)
            self.assertGreater(sim.brinkman_hydrodynamic_factor, 0.0)
            self.assertLessEqual(sim.brinkman_hydrodynamic_factor, 1.0)
            self.assertLessEqual(sim.macroscopic_diffusivity_um2_s, sim.effective_pore_diffusivity_um2_s)
            self.assertLessEqual(sim.effective_pore_diffusivity_um2_s, sim.stokes_free_diffusivity_um2_s)

    def test_05_zwitterionic_stealth_lnp_dornase_synergy(self) -> None:
        by_name = {p.name: simulate_mucus_filtration(p) for p in PROFILES}
        stealth = by_name["zwitterionic_peg_stealth_lnp_dornase"]
        dense_peg = by_name["dense_dmg_peg2000_lnp_cf_sputum"]

        self.assertGreater(stealth.macroscopic_diffusivity_um2_s, 4.0)
        self.assertGreater(stealth.epithelial_penetration_prob_pct, 80.0)
        self.assertGreater(stealth.ogston_steric_factor, dense_peg.ogston_steric_factor)
        self.assertGreater(stealth.brinkman_hydrodynamic_factor, dense_peg.brinkman_hydrodynamic_factor)

    def test_06_aav_capsid_vs_lnp_size_and_adhesion_tradeoff(self) -> None:
        by_name = {p.name: simulate_mucus_filtration(p) for p in PROFILES}
        aav = by_name["aav_capsid_25nm_cf_mucus"]
        stealth = by_name["zwitterionic_peg_stealth_lnp_dornase"]
        dense_peg = by_name["dense_dmg_peg2000_lnp_cf_sputum"]

        # 25 nm AAV has highest free and short-time pore diffusivity
        self.assertGreater(aav.stokes_free_diffusivity_um2_s, 20.0)
        self.assertGreater(aav.effective_pore_diffusivity_um2_s, stealth.effective_pore_diffusivity_um2_s)
        # Stronger capsid-mucin adhesion reduces mobile fraction below stealth LNP
        self.assertLess(aav.kramers_mobile_fraction, stealth.kramers_mobile_fraction)
        self.assertGreater(aav.epithelial_penetration_prob_pct, dense_peg.epithelial_penetration_prob_pct)

    def test_07_dense_peg_vs_low_peg_and_cationic_lipoplex(self) -> None:
        by_name = {p.name: simulate_mucus_filtration(p) for p in PROFILES}
        dense_peg = by_name["dense_dmg_peg2000_lnp_cf_sputum"]
        low_peg = by_name["low_peg_lnp_hyperconcentrated_cf"]
        cationic = by_name["unshielded_cationic_lipoplex_cf"]

        self.assertGreater(dense_peg.epithelial_penetration_prob_pct, 30.0)
        self.assertLess(low_peg.epithelial_penetration_prob_pct, 1.0)
        self.assertGreater(low_peg.epithelial_penetration_prob_pct, cationic.epithelial_penetration_prob_pct)
        self.assertLess(cationic.epithelial_penetration_prob_pct, 0.01)

    def test_08_poisoned_super_stokes_negative_control(self) -> None:
        poisoned = MucusFiltrationProfile(
            name="poisoned_super_stokes",
            label="Poisoned Super-Stokes",
            vector_class="Faulty Vector",
            mucus_condition="Negative Permeability",
            particle_radius_nm=20.0,
            fiber_radius_nm=2.0,
            mesh_pore_size_nm=400.0,
            brinkman_permeability_nm2=-500.0,
            surface_zeta_potential_mv=0.0,
            adhesive_barrier_kbt=-3.0,
            k_attempt_s1=120.0,
            k_trap_capture_s1=5.0,
            mucus_thickness_um=35.0,
            mcc_clearance_window_s=900.0,
        )
        sim = simulate_mucus_filtration(poisoned)
        self.assertFalse(sim.is_stokes_brinkman_bounded)
        self.assertFalse(sim.is_kramers_barrier_passive)
        self.assertGreater(sim.hindered_diffusivity_ratio_deff_d0, 1.0)

    def test_09_poisoned_cationic_trap_immobilization(self) -> None:
        poisoned_trap = MucusFiltrationProfile(
            name="poisoned_cationic_trap",
            label="Poisoned Cationic Trap",
            vector_class="High-Zeta Polyplex",
            mucus_condition="Dense CF Mucus",
            particle_radius_nm=75.0,
            fiber_radius_nm=4.5,
            mesh_pore_size_nm=90.0,
            brinkman_permeability_nm2=300.0,
            surface_zeta_potential_mv=45.0,
            adhesive_barrier_kbt=12.0,
            k_attempt_s1=120.0,
            k_trap_capture_s1=120.0,
            mucus_thickness_um=60.0,
            mcc_clearance_window_s=900.0,
        )
        sim = simulate_mucus_filtration(poisoned_trap)
        self.assertLess(sim.kramers_mobile_fraction, 1e-5)
        self.assertLess(sim.epithelial_penetration_prob_pct, 1e-10)

    def test_10_monotonic_mucus_thickness_and_barrier_scaling(self) -> None:
        base = PROFILES[0]
        sim_thin = simulate_mucus_filtration(base)
        thicker = MucusFiltrationProfile(
            **{**base.__dict__, "name": "thicker_mucus", "mucus_thickness_um": base.mucus_thickness_um * 2.0}
        )
        sim_thick = simulate_mucus_filtration(thicker)
        self.assertAlmostEqual(sim_thick.diffusive_transit_time_s, sim_thin.diffusive_transit_time_s * 4.0, places=5)
        self.assertLess(sim_thick.epithelial_penetration_prob_pct, sim_thin.epithelial_penetration_prob_pct)


if __name__ == "__main__":
    unittest.main()
