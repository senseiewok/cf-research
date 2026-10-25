#!/usr/bin/env python3
"""Unit tests for AlphaEvolve-CassetteShield (T-0433)."""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from alpha_evolve_cassette_shield import (
    UNSHIELDED_SUPEREXON_HOTSPOT_CDS,
    evaluate_cassette_shield,
    evolve_shielded_cds,
    translate_cds,
)


class TestAlphaEvolveCassetteShield(unittest.TestCase):
    def test_unshielded_cassette_exhibits_cpg_g4_and_microhomology_hazards(self) -> None:
        """Verify the baseline unshielded sequence contains the three AAV truncation/silencing hazards."""
        comp = evaluate_cassette_shield()
        self.assertGreaterEqual(comp.unshielded_cpg_count, 10)
        self.assertGreaterEqual(comp.unshielded_max_g4hunter, 1.0)
        self.assertGreaterEqual(comp.unshielded_max_repeat_bp, 10)

    def test_shielded_cassette_passes_all_six_formal_invariants(self) -> None:
        """Verify evolved sequence achieves 100% AA identity, zero CpG, safe G4, <=7 bp repeats, and tAI pause conservation."""
        comp = evaluate_cassette_shield()
        self.assertTrue(comp.all_invariants_passed, f"Unexpected unsat_core: {comp.unsat_core}")
        self.assertEqual(comp.amino_acid_identity_pct, 100.0)
        self.assertEqual(comp.shielded_cpg_count, 0)
        self.assertLess(comp.shielded_max_g4hunter, 0.85)
        self.assertLessEqual(comp.shielded_max_repeat_bp, 7)
        self.assertTrue(comp.pause_valley_preserved)
        self.assertLessEqual(comp.max_tai_profile_drift, 0.12)
        self.assertGreater(
            comp.estimated_intact_aav_genome_yield_shielded_pct,
            comp.estimated_intact_aav_genome_yield_unshielded_pct + 25.0,
        )

    def test_naive_max_codon_optimization_triggers_tai_pause_violation(self) -> None:
        """Replacing the inter-domain pause valley with rapid CTG/ATC codons triggers UNSAT_TAI_PAUSE_DRIFT."""
        shielded = evolve_shielded_cds(UNSHIELDED_SUPEREXON_HOTSPOT_CDS)
        # Corrupt the pause valley (codons 33..44, nt 99..132) with fast synonymous codons for L,L,I,V,S,C,T,N,H,D,I
        fast_valley = "CTGCTGATCGTGTCCTGCACCAACCACGACATC"
        corrupted = shielded[:99] + fast_valley + shielded[132:]
        self.assertEqual(translate_cds(corrupted), translate_cds(UNSHIELDED_SUPEREXON_HOTSPOT_CDS))
        comp_bad = evaluate_cassette_shield(UNSHIELDED_SUPEREXON_HOTSPOT_CDS, corrupted)
        self.assertFalse(comp_bad.inv_tai_profile_passed)
        self.assertFalse(comp_bad.all_invariants_passed)
        self.assertTrue(any("UNSAT_TAI_PAUSE_DRIFT" in u for u in comp_bad.unsat_core))

    def test_short_cds_and_rna_alphabet_normalization_succeed(self) -> None:
        """Verify RNA ('U') transcripts and short (<35-codon) sequences do not raise KeyError or empty max() ValueError."""
        rna_short = "AUGUGCCGCGGCGGCCGCGUGGGCGAGCGCCUG"  # 11 codons with U, CpGs, and G-runs
        opt_short = evolve_shielded_cds(rna_short)
        self.assertEqual(translate_cds(opt_short), translate_cds(rna_short))
        comp_short = evaluate_cassette_shield(rna_short.replace("U", "T"), opt_short)
        self.assertEqual(comp_short.amino_acid_identity_pct, 100.0)
        self.assertEqual(comp_short.shielded_cpg_count, 0)


if __name__ == "__main__":
    unittest.main()
