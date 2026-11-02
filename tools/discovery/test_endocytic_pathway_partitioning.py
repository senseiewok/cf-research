"""Unit tests for Macropinocytosis vs Clathrin-Mediated Endocytosis Flux in Stem Cells (T-0561).

Primary Biophysical Grounding & Evidence:
  - DOI:10.1038/s41565-020-0669-6 (Selective organ targeting SORT lipid nanoparticles for tissue-specific mRNA delivery)
  - DOI:10.1038/s41587-023-01861-1 (Inhaled delivery of lipid nanoparticles to respiratory epithelia)
  - DOI:10.1152/physrev.00004.2021 (Physiology and pathophysiology of human airway mucus)
  - NCT05668741 (Phase 1/2 VX-522 inhaled mRNA LNP study in CF)
"""
from __future__ import annotations

import unittest
from tools.discovery.endocytic_pathway_partitioning import (
    EndocyticPathwayParameters,
    compute_clathrin_mediated_flux,
    compute_macropinocytic_flux,
    compute_recycling_fraction,
    evaluate_endocytic_pathway_partitioning,
    falsify_endocytic_pathway_partitioning,
    simulate_endocytic_partitioning,
)


class EndocyticPathwayPartitioningTests(unittest.TestCase):
    """Test suite verifying CME and macropinocytic uptake kinetics, endosomal sorting, and stem cell transfection."""

    def setUp(self) -> None:
        self.report = evaluate_endocytic_pathway_partitioning()

    def test_report_structure_and_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0561")
        self.assertTrue(self.report["endocytic_partitioning_verified"])
        self.assertIn("Basal_Stem_SORT_Optimized", self.report["profiles"])
        self.assertIn("Basal_Stem_Standard_LNP", self.report["profiles"])
        self.assertIn("Quiescent_Stem_Untargeted", self.report["profiles"])

    def test_clathrin_mediated_flux_saturation(self) -> None:
        # J_CME = V_max * C / (Km + C)
        # At C = Km = 25 nM, V_max = 120: J_CME = 120 * 25 / 50 = 60.0 particles/cell/min
        j_half = compute_clathrin_mediated_flux(25.0, v_max=120.0, km_nm=25.0)
        self.assertAlmostEqual(j_half, 60.0, places=2)

        # At C = 0: J_CME = 0
        self.assertEqual(compute_clathrin_mediated_flux(0.0), 0.0)

        # High concentration approaches V_max (e.g. C = 1000 nM -> 120 * 1000 / 1025 = 117.07)
        j_high = compute_clathrin_mediated_flux(1000.0, v_max=120.0, km_nm=25.0)
        self.assertGreater(j_high, 115.0)
        self.assertLessEqual(j_high, 120.0)

    def test_macropinocytic_flux_scaling(self) -> None:
        # J_macro = k_macro * C * Phi_ruffle
        # At C = 30 nM, k_macro = 3.5, Phi_ruffle = 0.85: J_macro = 3.5 * 30 * 0.85 = 89.25
        j_macro = compute_macropinocytic_flux(30.0, k_macro=3.5, ruffling_index=0.85)
        self.assertAlmostEqual(j_macro, 89.25, places=2)

        # Zero ruffling yields zero macropinocytosis
        self.assertEqual(compute_macropinocytic_flux(30.0, ruffling_index=0.0), 0.0)

    def test_endosomal_sorting_and_recycling(self) -> None:
        # f_recycle = k_recycle / (k_sort + k_recycle)
        # For k_sort = 0.040, k_recycle = 0.020: total = 0.060, f_recycle = 0.020 / 0.060 = 0.3333
        f_rec = compute_recycling_fraction(0.040, 0.020)
        self.assertAlmostEqual(f_rec, 1.0 / 3.0, places=4)

    def test_transfection_yield_stratification(self) -> None:
        sort_opt = self.report["profiles"]["Basal_Stem_SORT_Optimized"]
        std_lnp = self.report["profiles"]["Basal_Stem_Standard_LNP"]
        quiescent = self.report["profiles"]["Quiescent_Stem_Untargeted"]
        cme_inhib = self.report["profiles"]["Clathrin_Inhibited_Rescue"]

        # SORT-optimized: robust stem cell transfection (>= 500 copies/cell)
        self.assertGreaterEqual(sort_opt["cytosolic_yield_mrna_per_cell"], 500.0)
        self.assertEqual(sort_opt["transfection_status"], "ROBUST_STEM_CELL_TRANSFECTION")

        # Standard LNP: moderate partial transfection (100 - 500 copies/cell)
        self.assertGreaterEqual(std_lnp["cytosolic_yield_mrna_per_cell"], 100.0)
        self.assertLess(std_lnp["cytosolic_yield_mrna_per_cell"], 500.0)
        self.assertEqual(std_lnp["transfection_status"], "MODERATE_PARTIAL_TRANSFECTION")

        # Quiescent untargeted: sub-therapeutic delivery failure (< 100 copies/cell)
        self.assertLess(quiescent["cytosolic_yield_mrna_per_cell"], 100.0)
        self.assertEqual(quiescent["transfection_status"], "SUB_THERAPEUTIC_DELIVERY_FAILURE")

        # Clathrin-inhibited: macropinocytosis compensates when actin ruffling is stimulated
        self.assertEqual(cme_inhib["j_cme_particles_per_min"], 0.0)
        self.assertEqual(cme_inhib["macro_flux_fraction"], 1.0)
        self.assertGreater(cme_inhib["cytosolic_yield_mrna_per_cell"], 300.0)

    def test_ah10_dynamic_falsification_negative_control(self) -> None:
        fals_res = falsify_endocytic_pathway_partitioning()
        self.assertEqual(fals_res["rule"], "AH-10")
        self.assertEqual(fals_res["verdict"], "PASS")
        self.assertTrue(fals_res["negative_control_passed"])
        self.assertEqual(fals_res["internalization_flux"], 0.0)
        self.assertEqual(fals_res["cytosolic_yield"], 0.0)


if __name__ == "__main__":
    unittest.main()
