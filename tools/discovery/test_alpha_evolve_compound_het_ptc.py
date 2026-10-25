#!/usr/bin/env python3
"""Unit tests for AlphaEvolve-CompoundHet-PTC (T-0432)."""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from alpha_evolve_compound_het_ptc import (
    ALLELE_BIOPHYSICAL_CATALOG,
    CompoundHetRegimen,
    evaluate_compound_het_genotype,
    optimize_compound_het_pareto_regimen,
)


class TestAlphaEvolveCompoundHetPtc(unittest.TestCase):
    def test_solo_readthrough_fails_15pct_threshold_due_to_neo_missense_misfolding(self) -> None:
        """Solo small-molecule readthrough (ELX-02-like) yields <3.5% WT current without ETI salvage."""
        res = evaluate_compound_het_genotype(
            "G542X",
            "W1282X",
            CompoundHetRegimen(
                small_mol_readthrough_u=0.80,
                ace_trna_cognate_u=0.0,
                nmd_stabilization_u=0.0,
                eti_corrector_u=0.0,
                ivacaftor_potentiator_u=0.0,
            ),
        )
        self.assertTrue(res.all_invariants_passed)
        self.assertLess(res.total_current_pct_wt, 5.0)
        self.assertFalse(res.crosses_15pct_responder_threshold)

    def test_readthrough_plus_eti_salvages_g542_and_w1282_neo_missense_proteoforms(self) -> None:
        """Co-administering ETI + NMD stabilization rescues G542W/C/R and W1282W/C/R above 15% WT."""
        summary = optimize_compound_het_pareto_regimen("G542X", "W1282X")
        self.assertLess(summary["solo_readthrough_pct_wt"], 5.0)
        self.assertGreaterEqual(summary["readthrough_plus_eti_pct_wt"], 15.0)
        self.assertGreater(
            summary["ace_trna_plus_eti_pct_wt"],
            summary["readthrough_plus_eti_pct_wt"],
        )
        self.assertGreaterEqual(summary["neo_missense_eti_epistatic_multiplier"], 4.0)

    def test_w1282x_Outperforms_g542x_under_near_cognate_uga_readthrough(self) -> None:
        """W1282X (native Trp) receives 52% native W1282W under UGA readthrough, exceeding G542X (0% native Gly)."""
        res = evaluate_compound_het_genotype(
            "G542X",
            "W1282X",
            CompoundHetRegimen(
                small_mol_readthrough_u=0.85,
                nmd_stabilization_u=0.65,
                eti_corrector_u=0.90,
                ivacaftor_potentiator_u=0.90,
            ),
        )
        g542_net = res.allele_1_breakdown["allele_net_current_pct_wt"]
        w1282_net = res.allele_2_breakdown["allele_net_current_pct_wt"]
        self.assertGreater(w1282_net, g542_net)
        self.assertEqual(res.allele_1_breakdown["proteoform_currents_pct_wt"]["G542G_WT"], 0.0)
        self.assertGreater(res.allele_2_breakdown["proteoform_currents_pct_wt"]["W1282W_WT"], 3.0)

    def test_compound_het_g542x_and_3849_splice_allele_synergy(self) -> None:
        """Compound heterozygote G542X / 3849+10kbC>T responds to combined ASO + readthrough + ETI."""
        res = evaluate_compound_het_genotype(
            "G542X",
            "3849+10kbC>T",
            CompoundHetRegimen(
                small_mol_readthrough_u=0.80,
                nmd_stabilization_u=0.50,
                eti_corrector_u=0.85,
                ivacaftor_potentiator_u=0.85,
                splice_aso_u=0.85,
            ),
        )
        self.assertTrue(res.all_invariants_passed)
        self.assertTrue(res.crosses_15pct_responder_threshold)
        self.assertEqual(len(res.ussing_trace_timepoints), 11)
        self.assertGreater(res.ussing_delta_isc_ua_cm2, 5.4)

    def test_all_catalog_alleles_conserve_probability_simplex_mass(self) -> None:
        """Verify INV-READTHROUGH-MASS across all catalog allele pairings."""
        allele_ids = sorted(ALLELE_BIOPHYSICAL_CATALOG.keys())
        for a1 in allele_ids:
            res = evaluate_compound_het_genotype(a1, "W1282X")
            self.assertTrue(res.inv_readthrough_mass_passed, f"Mass conservation failed for {a1}")
            self.assertTrue(res.inv_chaperone_bound_passed, f"Chaperone bound failed for {a1}")
            self.assertTrue(res.inv_biallelic_monotonicity_passed, f"Monotonicity failed for {a1}")

    def test_deconvoluted_proteoform_currents_sum_to_scaled_allele_net_current(self) -> None:
        """Verify ERAD-scaled proteoform currents + truncated current sum to allele_net_current_pct_wt."""
        res = evaluate_compound_het_genotype(
            "G542X",
            "W1282X",
            CompoundHetRegimen(
                small_mol_readthrough_u=0.90,
                nmd_stabilization_u=0.80,
                eti_corrector_u=0.90,
                ivacaftor_potentiator_u=0.90,
            ),
        )
        for breakdown in (res.allele_1_breakdown, res.allele_2_breakdown):
            sum_parts = (
                sum(breakdown["proteoform_currents_pct_wt"].values())
                + breakdown["truncated_current_pct_wt"]
            )
            self.assertAlmostEqual(
                sum_parts,
                breakdown["allele_net_current_pct_wt"],
                delta=0.015,
            )


if __name__ == "__main__":
    unittest.main()
