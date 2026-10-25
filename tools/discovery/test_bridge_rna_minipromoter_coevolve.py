#!/usr/bin/env python3
"""Unit tests for Pan-Mutation IS110 Bridge-RNA Recombinase & Airway Minipromoter Co-Evolution (T-0427)."""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from bridge_rna_minipromoter_coevolve import (
    BridgeRnaAavCassette,
    evaluate_bridge_rna_cassette,
    run_halton_patient_cohort_uq,
)


class TestBridgeRnaMinipromoterCoevolve(unittest.TestCase):
    def test_single_aav_payload_budget_invariant(self) -> None:
        ev = evaluate_bridge_rna_cassette()
        self.assertEqual(ev.total_payload_bp, 3630)
        self.assertEqual(ev.packaging_headroom_bp, 1070)
        self.assertTrue(ev.inv_aav_cap_passed)
        self.assertTrue(ev.all_invariants_passed)
        self.assertEqual(ev.unsat_core, [])

    def test_oversized_cassette_triggers_unsat_core(self) -> None:
        bad = BridgeRnaAavCassette(cftr_super_exon_15_27_bp=3400)
        ev = evaluate_bridge_rna_cassette(bad)
        self.assertFalse(ev.inv_aav_cap_passed)
        self.assertFalse(ev.all_invariants_passed)
        self.assertTrue(any("UNSAT_AAV_CAP" in s for s in ev.unsat_core))

    def test_monte_carlo_cohort_credible_bands_n256(self) -> None:
        uq = run_halton_patient_cohort_uq(256)
        self.assertGreaterEqual(uq["isc_p10"], 14.0)
        self.assertGreaterEqual(uq["isc_p50"], 18.0)
        self.assertGreaterEqual(uq["asl_p10"], 7.00)
        self.assertLessEqual(uq["lci_p50"], -1.85)


if __name__ == "__main__":
    unittest.main()
