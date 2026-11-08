"""Unit tests for Round 117: Hilbert's Sixth Problem Boltzmann-to-Navier-Stokes Limit & Mucus Knudsen Slip Bridge."""

import unittest
from tools.discovery.round117_hilbert_sixth_boltzmann_navier_stokes_mucus_bridge import (
    CITATIONS,
    HYPOTHESIS,
    PARAMETER_PROVENANCE,
    PRIMARY_DOI,
    PROFILES,
    ROUND_ID,
    SECONDARY_DOIS,
    TASK_ID,
    TITLE,
    HilbertSixthMucusProfile,
    compute_cunningham_knudsen_slip,
    run_falsification,
    run_falsification_checks,
    run_round117_verification,
    simulate_hilbert_sixth_mucus_bridge,
)


class TestRound117HilbertSixthBoltzmannNavierStokesMucusBridge(unittest.TestCase):
    def test_01_metadata_and_ah10_contract(self) -> None:
        self.assertEqual(ROUND_ID, "round117")
        self.assertEqual(TASK_ID, "T-0753")
        self.assertEqual(PRIMARY_DOI, "10.1007/s00222-024-01298-w")
        self.assertEqual(
            SECONDARY_DOIS,
            ["10.1126/science.1223012", "10.1016/j.addr.2008.11.002"],
        )
        self.assertIn("Hilbert's Sixth Problem", TITLE)
        self.assertIn("Knudsen", TITLE)
        self.assertIn("Deng-Hani-Ma", HYPOTHESIS)
        self.assertEqual(len(CITATIONS), 3)
        self.assertTrue(all(c.get("scite_verified") is True for c in CITATIONS))

        ah10 = run_falsification()
        self.assertTrue(ah10["passed"])
        self.assertEqual(ah10["rule"], "AH-10")
        self.assertEqual(ah10["round_id"], ROUND_ID)
        self.assertEqual(ah10["task_id"], TASK_ID)

    def test_02_parameter_provenance_bounds(self) -> None:
        for p_name, p_meta in PARAMETER_PROVENANCE.items():
            self.assertIn("unit", p_meta, msg=f"Missing unit for {p_name}")
            self.assertIn("kind", p_meta, msg=f"Missing kind for {p_name}")
            self.assertIn("source", p_meta, msg=f"Missing source for {p_name}")
            self.assertIn("bounds", p_meta, msg=f"Missing bounds for {p_name}")
            low, high = p_meta["bounds"]
            self.assertLess(low, high)
            if "value" in p_meta:
                val = p_meta["value"]
                self.assertGreaterEqual(val, low)
                self.assertLessEqual(val, high)

    def test_03_five_mucus_profiles_inventory(self) -> None:
        expected_names = [
            "healthy_pcl_peg_lnp_boltzmann_slip",
            "f508del_eti_rescued_peg_lnp",
            "class1_zwitterionic_mpp_hypertonic_saline",
            "class1_dehydrated_cf_sputum_standard_lnp",
            "severe_cf_plug_unshielded_cationic_lnp",
        ]
        self.assertEqual([p.name for p in PROFILES], expected_names)

    def test_04_all_falsification_checks_and_verification_contract(self) -> None:
        report = run_falsification_checks()
        self.assertTrue(report["passed"])
        for check_name, check_val in report["checks"].items():
            self.assertTrue(check_val, msg=f"Check failed: {check_name}")

        ver = run_round117_verification()
        self.assertEqual(ver["round"], 117)
        self.assertEqual(ver["status"], "PASS")
        self.assertEqual(ver["doi"], "10.1007/s00222-024-01298-w")
        self.assertEqual(len(ver["secondary_dois"]), 2)
        self.assertEqual(len(ver["genotypes"]), 5)
        self.assertEqual(len(ver["negative_controls"]), 3)
        self.assertTrue(all(nc["rejected"] for nc in ver["negative_controls"]))
        self.assertTrue(ver["summary_metrics"]["parametric_sensitivity_25pct_passed"])

    def test_05_cunningham_knudsen_slip_and_osmotic_brush_invariants(self) -> None:
        self.assertEqual(compute_cunningham_knudsen_slip(0.0), 1.0)
        for p in PROFILES:
            res = simulate_hilbert_sixth_mucus_bridge(p)
            self.assertGreater(res.pore_knudsen_number, 0.0)
            self.assertGreaterEqual(res.cunningham_knudsen_slip_factor, 1.0)
            self.assertGreater(res.hilbert6_continuum_defect_norm, 0.0)
            self.assertLess(res.hilbert6_continuum_defect_norm, 1.0)
            self.assertGreater(res.pcl_brush_height_um, 0.0)
            self.assertLessEqual(res.pcl_brush_height_um, 7.0)
            self.assertTrue(res.is_cunningham_slip_valid)
            self.assertTrue(res.is_diffusivity_bound_valid)
            self.assertTrue(res.is_osmotic_pcl_bound_valid)

    def test_06_monotonic_ordering_and_negative_controls(self) -> None:
        results = {p.name: simulate_hilbert_sixth_mucus_bridge(p) for p in PROFILES}
        h = results["healthy_pcl_peg_lnp_boltzmann_slip"]
        eti = results["f508del_eti_rescued_peg_lnp"]
        mpp = results["class1_zwitterionic_mpp_hypertonic_saline"]
        dehyd = results["class1_dehydrated_cf_sputum_standard_lnp"]
        plug = results["severe_cf_plug_unshielded_cationic_lnp"]

        self.assertGreater(h.effective_diffusivity_ratio, eti.effective_diffusivity_ratio)
        self.assertGreater(eti.effective_diffusivity_ratio, mpp.effective_diffusivity_ratio)
        self.assertGreater(mpp.effective_diffusivity_ratio, dehyd.effective_diffusivity_ratio)
        self.assertGreater(dehyd.effective_diffusivity_ratio, plug.effective_diffusivity_ratio)

        self.assertGreater(h.translocation_probability, eti.translocation_probability)
        self.assertGreater(eti.translocation_probability, mpp.translocation_probability)
        self.assertGreater(mpp.translocation_probability, dehyd.translocation_probability)
        self.assertGreater(dehyd.translocation_probability, plug.translocation_probability)

        bad_slip = simulate_hilbert_sixth_mucus_bridge(
            HilbertSixthMucusProfile(
                name="bad_slip",
                description="bad",
                mucin_solids_wt_pct=4.8,
                mesh_correlation_nm=125.0,
                lnp_radius_nm=35.0,
                zeta_abs_mv=1.8,
                interstitial_eta_rel=11.0,
                poison_mode="sub_unitary_cunningham_slip_violation",
            )
        )
        self.assertFalse(bad_slip.is_cunningham_slip_valid)


if __name__ == "__main__":
    unittest.main()
