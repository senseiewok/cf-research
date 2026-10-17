#!/usr/bin/env python3
"""Unit tests for cross_domain_miner.py (T-0179)."""
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
            macrophage_stealth_score=0.85
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


if __name__ == "__main__":
    unittest.main()
