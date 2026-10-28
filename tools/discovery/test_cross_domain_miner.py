#!/usr/bin/env python3
"""Unit tests for cross_domain_miner.py and CrossDiseaseMechanisticBridge (T-0179)."""
import unittest
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import cross_domain_miner as cdm


class TestCrossDomainMiner(unittest.TestCase):
    def test_curated_innovations_integrity(self):
        self.assertGreaterEqual(len(cdm.CURATED_INNOVATIONS), 5)
        for item in cdm.CURATED_INNOVATIONS:
            self.assertTrue(item["id"].startswith("CDI-"))
            self.assertTrue(item["primary_doi"])
            self.assertTrue(item["donor_field"])
            self.assertGreater(item["particle_diameter_nm"], 0.0)
            self.assertLessEqual(abs(item["zeta_potential_mv"]), 50.0)

    def test_score_steric_mesh_penetration(self):
        # Optimal size 50 nm should score high
        self.assertGreaterEqual(cdm.score_steric_mesh_penetration(50.0), 0.90)
        # Oversized particle > 140 nm should score low
        self.assertLess(cdm.score_steric_mesh_penetration(150.0), 0.20)
        # Non-positive diameter should return 0.0
        self.assertEqual(cdm.score_steric_mesh_penetration(0.0), 0.0)

    def test_score_zeta_potential_neutrality(self):
        # Near neutral (-1 mV) should score near 1.0
        self.assertGreaterEqual(cdm.score_zeta_potential_neutrality(-1.0), 0.95)
        # Highly positive (+30 mV) should receive heavy penalty
        self.assertLess(cdm.score_zeta_potential_neutrality(30.0), 0.30)

    def test_calculate_bfi(self):
        bfi = cdm.calculate_biophysical_feasibility_index(
            diameter_nm=60.0,
            zeta_mv=-1.5,
            endosomal_escape_score=0.90,
            macrophage_stealth_score=0.85,
        )
        self.assertGreaterEqual(bfi, 0.85)
        self.assertLessEqual(bfi, 1.0)

    def test_screen_innovations_filtering(self):
        # All innovations should meet min_bfi 0.70
        all_results = cdm.screen_innovations(min_bfi=0.70)
        self.assertEqual(len(all_results), len(cdm.CURATED_INNOVATIONS))

        # Filter by specific field
        bbb_results = cdm.screen_innovations(donor_field="blood_brain_barrier")
        self.assertEqual(len(bbb_results), 1)
        self.assertEqual(bbb_results[0]["id"], "CDI-002")

        # Ranking is descending by BFI
        bfi_list = [r["bfi"] for r in all_results]
        self.assertEqual(bfi_list, sorted(bfi_list, reverse=True))

    def test_cross_disease_bridges_invariants(self):
        audited = cdm.CrossDiseaseMechanisticBridge.audit_all_curated_bridges()
        self.assertEqual(len(audited), 4)
        lead_ids = [b["lead_id"] for b in audited]
        self.assertEqual(lead_ids, ["L-0039", "L-0040", "L-0041", "L-0042"])
        for b in audited:
            self.assertTrue(b["verification"]["inv_cf_anchor_satisfied"])
            self.assertTrue(b["verification"]["inv_mech_isomorphism_satisfied"])
            self.assertFalse(b["verification"]["clinical_extrapolation_permitted"])

        # Negative control: clinical_extrapolation_permitted=True must raise ValueError
        bad_bridge = dict(cdm.CURATED_CROSS_DISEASE_BRIDGES[0])
        bad_bridge["clinical_extrapolation_permitted"] = True
        with self.assertRaises(ValueError):
            cdm.CrossDiseaseMechanisticBridge.verify_bridge_invariants(bad_bridge)

        # Negative control: unknown cf_bottleneck must raise ValueError
        bad_anchor = dict(cdm.CURATED_CROSS_DISEASE_BRIDGES[0])
        bad_anchor["cf_bottleneck"] = "unanchored_non_cf_topic"
        with self.assertRaises(ValueError):
            cdm.CrossDiseaseMechanisticBridge.verify_bridge_invariants(bad_anchor)

    def test_abc_transporter_folding_and_negative_control(self):
        # Valid domain-interface folding rescue
        res_valid = cdm.CrossDiseaseMechanisticBridge.evaluate_abc_transporter_folding_bridge(
            delta_delta_g_interface_kcal_mol=3.2,
            corrector_stabilization_kcal_mol=-2.6,
            rnf5_erad_escape_factor=1.4,
            claims_channel_gating_equivalence=False,
        )
        self.assertTrue(res_valid["isomorphic"])
        self.assertGreater(res_valid["escaped_folded_fraction"], 0.30)

        # Negative control: equating CFTR channel Po to active pump turnover k_cat must fail
        res_invalid = cdm.CrossDiseaseMechanisticBridge.evaluate_abc_transporter_folding_bridge(
            delta_delta_g_interface_kcal_mol=3.2,
            corrector_stabilization_kcal_mol=-2.6,
            claims_channel_gating_equivalence=True,
        )
        self.assertFalse(res_invalid["isomorphic"])
        self.assertEqual(res_invalid["escaped_folded_fraction"], 0.0)
        self.assertIn("REJECTED", res_invalid["falsification_reason"])

    def test_ptc_readthrough_folding_salvage(self):
        # Standalone readthrough with severe near-cognate cavity penalty (+2.8 kcal/mol)
        standalone = cdm.CrossDiseaseMechanisticBridge.evaluate_ptc_readthrough_folding_salvage(
            nmd_escape_fraction=0.65,
            basal_readthrough_prob=0.12,
            near_cognate_cavity_penalty_kcal_mol=2.8,
            corrector_rescue_kcal_mol=0.0,
        )
        self.assertTrue(standalone["trimodal_salvage_required"])

        # Trimodal salvage with corrector rescue (-2.5 kcal/mol) increases functional salvage >2.5x
        trimodal = cdm.CrossDiseaseMechanisticBridge.evaluate_ptc_readthrough_folding_salvage(
            nmd_escape_fraction=0.65,
            basal_readthrough_prob=0.12,
            near_cognate_cavity_penalty_kcal_mol=2.8,
            corrector_rescue_kcal_mol=-2.5,
        )
        self.assertGreater(
            trimodal["functional_salvage_fraction"],
            standalone["functional_salvage_fraction"] * 2.5,
        )

    def test_cotranslational_pause_bridge_negative_control(self):
        # Preserved pause valleys (r = 0.92, dwell = 1.6 s)
        preserved = cdm.CrossDiseaseMechanisticBridge.evaluate_cotranslational_pause_bridge(
            pause_valley_pearson_r=0.92,
            interdomain_dwell_sec=1.6,
            global_cai_maximized=False,
        )
        self.assertTrue(preserved["pause_valley_preserved"])
        self.assertEqual(preserved["inv_tai_profile_status"], "PASS")

        # Negative control: naive global CAI maximization (CAI = 1.0) triggers collision risk
        naive_cai = cdm.CrossDiseaseMechanisticBridge.evaluate_cotranslational_pause_bridge(
            pause_valley_pearson_r=0.92,
            interdomain_dwell_sec=1.6,
            global_cai_maximized=True,
        )
        self.assertFalse(naive_cai["pause_valley_preserved"])
        self.assertEqual(naive_cai["inv_tai_profile_status"], "COLLISION_RISK")

    def test_muco_obstructive_rheology_negative_control(self):
        # Negative control: rhDNase in low-DNA non-CF plug (0.2 mg/mL DNA) must be rejected
        low_dna = cdm.CrossDiseaseMechanisticBridge.evaluate_muco_obstructive_rheology_bridge(
            solids_wt_pct=6.5,
            extracellular_dna_mg_ml=0.2,
            modality="rhDNase",
        )
        self.assertFalse(low_dna["mechanistic_transfer_valid"])
        self.assertIn("REJECTED", low_dna["falsification_reason"])

        # Purulent high-DNA CF plug (3.2 mg/mL DNA) with rhDNase succeeds
        high_dna = cdm.CrossDiseaseMechanisticBridge.evaluate_muco_obstructive_rheology_bridge(
            solids_wt_pct=6.5,
            extracellular_dna_mg_ml=3.2,
            modality="rhDNase",
        )
        self.assertTrue(high_dna["mechanistic_transfer_valid"])


if __name__ == "__main__":
    unittest.main()
