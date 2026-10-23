#!/usr/bin/env python3
"""Unit tests for CFTR Cryptic Splice Site Pseudoexon U1 snRNP Steric Decoy Thermodynamics Model (T-0352).

Acceptance Criteria:
  python tools/discovery/u1_snrnp_decoy_thermodynamics.py --check
  pytest tools/discovery/test_u1_snrnp_decoy_thermodynamics.py
"""

from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from u1_snrnp_decoy_thermodynamics import (
    DECOYS_CATALOG,
    MUTATIONS_CATALOG,
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
    compute_competitive_occupancy,
    compute_nearest_neighbor_thermodynamics,
    run_checks,
    simulate_titration_curve,
)


class TestU1snRNPDecoyThermodynamics(unittest.TestCase):
    """Unit test suite for T-0352."""

    def test_catalog_entries(self):
        self.assertIn("3849_plus_10kb_C_to_T", MUTATIONS_CATALOG)
        self.assertIn("1811_plus_1_6kb_A_to_G", MUTATIONS_CATALOG)
        self.assertIn("SPL84_2MOE", DECOYS_CATALOG)
        self.assertIn("PMO_StericBlocker", DECOYS_CATALOG)
        self.assertIn("ExSpeU1_Engineered", DECOYS_CATALOG)

    def test_regulatory_disclaimers(self):
        self.assertIn("ClinicalTrials.gov does not evaluate safety or efficacy", NLM_ATTRIBUTION)
        self.assertIn("Preclinical in-silico biophysical modeling", RESEARCH_BOUNDARY)
        for prohibited in ["proven to cure", "guarantees clinical benefit", "cures cystic fibrosis"]:
            self.assertNotIn(prohibited, NLM_ATTRIBUTION.lower())
            self.assertNotIn(prohibited, RESEARCH_BOUNDARY.lower())

    def test_nearest_neighbor_thermodynamics(self):
        wt = compute_nearest_neighbor_thermodynamics("GCAAGU")
        mut = compute_nearest_neighbor_thermodynamics("GUAAGU")

        self.assertAlmostEqual(wt.delta_g_37_kcal_mol, -3.2, places=1)
        self.assertAlmostEqual(mut.delta_g_37_kcal_mol, -10.4, places=1)
        self.assertLess(mut.delta_g_37_kcal_mol, wt.delta_g_37_kcal_mol)
        self.assertLess(mut.kd_nm, wt.kd_nm)
        self.assertLess(mut.kd_nm, 100.0)  # Nanomolar affinity for consensus
        self.assertGreater(wt.kd_nm, 10000.0)  # Micromolar affinity for WT non-consensus

    def test_competitive_mass_action(self):
        # At zero dose, U1 occupies majority of donor sites
        t_d0, t_u1_0 = compute_competitive_occupancy(0.0, 0.75, 30.0, 25.0)
        self.assertAlmostEqual(t_d0, 0.0, places=4)
        self.assertGreater(t_u1_0, 0.40)

        # At high dose (100 nM), decoy displaces U1
        t_d100, t_u1_100 = compute_competitive_occupancy(100.0, 0.75, 30.0, 25.0)
        self.assertGreater(t_d100, 0.85)
        self.assertLess(t_u1_100, 0.10)
        self.assertLessEqual(t_d100 + t_u1_100, 1.0)

    def test_dose_response_and_genotype_contrast(self):
        mut = MUTATIONS_CATALOG["3849_plus_10kb_C_to_T"]
        decoy = DECOYS_CATALOG["SPL84_2MOE"]
        titration = simulate_titration_curve(mut, decoy)

        self.assertEqual(len(titration), 10)
        p0 = titration[0]
        p_hi = titration[-1]

        self.assertEqual(p0.decoy_dose_nm, 0.0)
        self.assertEqual(p0.pseudoexon_exclusion_pct, 0.0)
        self.assertGreater(p_hi.pseudoexon_exclusion_pct, 85.0)
        self.assertGreater(p_hi.full_length_mrna_pct_homozygous, p_hi.full_length_mrna_pct_compound_het)
        self.assertTrue(p_hi.curative_benchmark_met)

    def test_run_checks_deterministic(self):
        code = run_checks()
        self.assertEqual(code, 0)


if __name__ == "__main__":
    unittest.main()
