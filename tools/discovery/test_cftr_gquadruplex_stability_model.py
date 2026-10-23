#!/usr/bin/env python3
"""Unit tests for CFTR mRNA Coding Sequence G-Quadruplex Stability & Ribosome Stalling Model (T-0316).

Acceptance Criteria:
  python -m unittest tools.discovery.test_cftr_gquadruplex_stability_model; python tools/discovery/cftr_gquadruplex_stability_model.py --check
"""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from cftr_gquadruplex_stability_model import (
    CFTR_DOMAINS,
    KNOWN_CFTR_RG4_MOTIFS,
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
    calculate_ribosome_elongation_kinetics,
    compute_g4hunter_score,
    get_domain_for_codon,
    run_checks,
    simulate_small_molecule_perturbation,
)


class TestCFTRGQuadruplexStabilityModel(unittest.TestCase):
    """Unit test suite for T-0316."""

    def test_domain_partition(self):
        self.assertEqual(len(CFTR_DOMAINS), 5)
        self.assertEqual(get_domain_for_codon(1), "MSD1")
        self.assertEqual(get_domain_for_codon(388), "MSD1")
        self.assertEqual(get_domain_for_codon(508), "NBD1")
        self.assertEqual(get_domain_for_codon(700), "R_DOMAIN")
        self.assertEqual(get_domain_for_codon(1000), "MSD2")
        self.assertEqual(get_domain_for_codon(1282), "NBD2")

    def test_g4hunter_propensity_scoring(self):
        # Canonical rG4 sequences should score well above 1.0
        g4_seq = "GGGAGGGUGGAGGCUUGGGG"
        non_g4_seq = "AUCGAUCGAUCGAUCG"
        self.assertGreater(compute_g4hunter_score(g4_seq), 1.5)
        self.assertLess(compute_g4hunter_score(non_g4_seq), 0.5)

    def test_known_cftr_motifs_thermodynamics(self):
        self.assertGreaterEqual(len(KNOWN_CFTR_RG4_MOTIFS), 5)
        for m in KNOWN_CFTR_RG4_MOTIFS:
            self.assertLess(m.delta_g_unfold_kcal_mol, -5.0)
            self.assertGreater(m.tm_c, 50.0)
            self.assertGreaterEqual(m.g4hunter_score, 1.0)

    def test_ribosome_dwell_kinetics_and_nbd1_pause(self):
        dwells = calculate_ribosome_elongation_kinetics()
        self.assertGreater(len(dwells), 0)

        # Confirm dwell times at rG4 are longer than baseline
        rg4_points = [p for p in dwells if p.rg4_present]
        self.assertGreater(len(rg4_points), 0)
        for p in rg4_points:
            self.assertGreater(p.perturbed_dwell_time_s, p.baseline_dwell_time_s)

        # NBD1 pre-F508 pause site validation
        nbd1_pauses = [p for p in rg4_points if p.domain == "NBD1" and p.codon_index < 508]
        self.assertTrue(len(nbd1_pauses) > 0)
        self.assertIn("productive_co_translational_folding", nbd1_pauses[0].co_translational_folding_role)

    def test_small_molecule_perturbation_and_dhx36_rescue(self):
        pert = simulate_small_molecule_perturbation()
        self.assertIn("Vehicle_Control", pert)
        self.assertIn("Pyridostatin_PDS", pert)
        self.assertIn("PDS_Plus_DHX36_Rescue", pert)

        veh = pert["Vehicle_Control"]
        pds = pert["Pyridostatin_PDS"]
        rescue = pert["PDS_Plus_DHX36_Rescue"]

        self.assertGreater(veh.cftr_protein_yield_pct, pds.cftr_protein_yield_pct)
        self.assertGreater(rescue.cftr_protein_yield_pct, pds.cftr_protein_yield_pct)
        self.assertAlmostEqual(veh.delta_tm_c, 0.0)
        self.assertGreater(pds.delta_tm_c, 10.0)

    def test_regulatory_disclaimers(self):
        self.assertIn("ClinicalTrials.gov does not evaluate safety or efficacy", NLM_ATTRIBUTION)
        self.assertIn("Preclinical in-silico biophysical modeling", RESEARCH_BOUNDARY)
        for prohibited in ["proven to cure", "guarantees clinical benefit", "cures cystic fibrosis"]:
            self.assertNotIn(prohibited, NLM_ATTRIBUTION.lower())
            self.assertNotIn(prohibited, RESEARCH_BOUNDARY.lower())

    def test_run_checks_acceptance(self):
        code = run_checks()
        self.assertEqual(code, 0)


if __name__ == "__main__":
    unittest.main()
