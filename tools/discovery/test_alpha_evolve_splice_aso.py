#!/usr/bin/env python3
"""Unit tests for AlphaEvolve-Splice (C2-v2, T-0434)."""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from alpha_evolve_splice_aso import (
    evaluate_aso_candidate,
    run_sliding_window_aso_walk,
)


class TestAlphaEvolveSpliceAso(unittest.TestCase):
    def test_sliding_window_walk_identifies_high_potency_invariant_safe_aso(self) -> None:
        """Verify the 1-bp walk across 3849+10kb C->T identifies candidates restoring >60% canonical WT mRNA."""
        walk = run_sliding_window_aso_walk()
        self.assertGreater(walk["total_windows_scanned"], 100)
        self.assertGreater(walk["invariant_verified_count"], 15)
        winner = walk["top_pareto_winner"]
        self.assertTrue(winner["all_invariants_passed"])
        self.assertGreaterEqual(winner["predicted_canonical_wt_mrna_pct"], 60.0)
        self.assertLessEqual(winner["dg_net_hybridization_kcal_mol"], -14.0)
        self.assertEqual(winner["cpg_count"], 0)
        self.assertEqual(winner["tlr78_gu_motif_count"], 0)
        self.assertFalse(winner["has_g_tetrad"])

    def test_distal_window_has_lower_splice_occlusion_than_cryptic_donor_center(self) -> None:
        """Windows far from the cryptic GU donor (offset 0) have lower occlusion than donor-spanning windows."""
        distal = evaluate_aso_candidate(start_offset_nt=0, length_nt=19)
        central = evaluate_aso_candidate(start_offset_nt=17, length_nt=19)
        self.assertGreater(
            central.donor_steric_occlusion_efficiency,
            distal.donor_steric_occlusion_efficiency + 0.40,
        )
        self.assertGreater(
            central.predicted_canonical_wt_mrna_pct,
            distal.predicted_canonical_wt_mrna_pct + 25.0,
        )

    def test_negative_control_immunogenic_and_g4_motifs_fail_invariants(self) -> None:
        """Verify that engineered CpG, TLR7/8 GUUG/UGUU motifs, and GGGG tetrads trigger formal UNSAT violations."""
        # Target containing CCCC (reverse complement GGGG) and CG + CAAC (reverse complement GUUG)
        bad_window = "UUUGAUUUGAUAAUGACCCCCGCAACUUGGUAAGUAACUAAUGCCAAUUGAAUAAAUUG"
        cand = evaluate_aso_candidate(start_offset_nt=12, length_nt=19, pre_mrna_window=bad_window)
        self.assertTrue(cand.has_g_tetrad)
        self.assertGreater(cand.cpg_count, 0)
        self.assertGreater(cand.tlr78_gu_motif_count, 0)
        self.assertFalse(cand.inv_no_g4_agg_passed)
        self.assertFalse(cand.inv_tlr_immune_safe_passed)
        self.assertFalse(cand.all_invariants_passed)
        self.assertTrue(any("UNSAT_G4_AGGREGATION" in u for u in cand.unsat_core))
        self.assertTrue(any("UNSAT_TLR_IMMUNE" in u for u in cand.unsat_core))

    def test_out_of_bounds_offset_raises_value_error(self) -> None:
        """Verify negative start offsets or slices exceeding pre_mrna_window length raise ValueError."""
        with self.assertRaises(ValueError):
            evaluate_aso_candidate(start_offset_nt=-2, length_nt=19)
        with self.assertRaises(ValueError):
            evaluate_aso_candidate(start_offset_nt=50, length_nt=19)


if __name__ == "__main__":
    unittest.main()
