#!/usr/bin/env python3
"""Unit tests for Airway Basal Stem Cell Tight Junction Kinetics (T-0429)."""
import math
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from airway_basal_tight_junction_kinetics import (
    CLINICAL_TRIAL_4D710,
    CLINICAL_TRIAL_VX522,
    NLM_ATTRIBUTION,
    PDB_ACTIVE_CFTR,
    PDB_DEPHOSPHO_CFTR,
    AirwayEpithelialLineageModel,
    TightJunctionKinetics,
    calculate_renkin_partitioning,
)


class TestAirwayBasalTightJunctionKinetics(unittest.TestCase):
    """Test suite verifying tight junction dynamics, Renkin partitioning, and clonal persistence."""

    def setUp(self) -> None:
        self.tj = TightJunctionKinetics(
            r0_nm=0.5,
            r_max_nm=35.0,
            teer0_ohm_cm2=500.0,
            tau_reseal_hours=2.0,
        )
        self.lineage = AirwayEpithelialLineageModel()

    def test_renkin_partitioning_boundaries(self) -> None:
        """Verify Renkin steric hindrance bounds at zero, intermediate, and cutoff ratios."""
        self.assertAlmostEqual(calculate_renkin_partitioning(0.0), 1.0, places=4)
        self.assertEqual(calculate_renkin_partitioning(1.0), 0.0)
        self.assertEqual(calculate_renkin_partitioning(1.2), 0.0)

        # Monotonic decrease in intermediate range
        p_25 = calculate_renkin_partitioning(0.25)
        p_50 = calculate_renkin_partitioning(0.50)
        p_75 = calculate_renkin_partitioning(0.75)
        self.assertGreater(p_25, p_50)
        self.assertGreater(p_50, p_75)
        self.assertGreater(p_75, 0.0)

    def test_tight_junction_teer_and_dilation_kinetics(self) -> None:
        """Verify TEER inverse-square relationship and exponential resealing kinetics."""
        self.assertAlmostEqual(self.tj.teer_at_radius(0.5), 500.0, places=2)
        # Dilated to 5.0 nm -> 100x drop in TEER
        self.assertAlmostEqual(self.tj.teer_at_radius(5.0), 5.0, places=2)

        # Dynamic pore opening during 2-hour duration
        r_during = self.tj.dynamic_pore_radius(1.0, duration_hours=2.0)
        self.assertEqual(r_during, 35.0)

        # Resealing at t = 2h + 2h (one half-life dt = 2h)
        r_reseal_1hl = self.tj.dynamic_pore_radius(4.0, duration_hours=2.0)
        expected_r = 0.5 + (35.0 - 0.5) * math.exp(-1.0)
        self.assertAlmostEqual(r_reseal_1hl, expected_r, places=2)

    def test_vector_paracellular_permeability(self) -> None:
        """Verify steric exclusion of AAV (r_h = 12.5 nm) at baseline and passage upon dilation."""
        # Baseline pore = 0.5 nm -> complete exclusion
        perm_baseline = self.tj.vector_paracellular_permeability(0.5, vector_radius_nm=12.5)
        self.assertEqual(perm_baseline, 0.0)

        # Dilated pore = 35.0 nm -> lambda = 12.5 / 35.0 = 0.357 -> permeable
        perm_dilated = self.tj.vector_paracellular_permeability(35.0, vector_radius_nm=12.5)
        self.assertGreater(perm_dilated, 0.10)
        self.assertLess(perm_dilated, 0.50)

    def test_apical_only_transduction_shedding_decay(self) -> None:
        """Verify that apical-only transduction clears within 180 days due to ciliated cell shedding."""
        sim = self.lineage.simulate_persistence(
            days=180.0,
            dt_days=1.0,
            initial_transduced_basal=0.0,
            initial_transduced_ciliated=250.0,
        )
        final_state = sim[-1]
        self.assertEqual(final_state["transduced_basal"], 0.0)
        self.assertEqual(final_state["transduced_ta"], 0.0)
        # Decays to less than 15% of initial by day 180 (exp(-180 / (60 / ln 2)) = 2^(-3) = 12.5%)
        self.assertLess(final_state["transduced_ciliated"], 35.0)

    def test_basal_transduction_clonal_durability(self) -> None:
        """Verify that transducing basal stem cells sustains ciliated expression at 365 days."""
        sim = self.lineage.simulate_persistence(
            days=365.0,
            dt_days=1.0,
            initial_transduced_basal=50.0,
            initial_transduced_ciliated=0.0,
        )
        final_state = sim[-1]
        # Basal stem cells self-renew, feeding differentiated pool
        self.assertGreater(final_state["transduced_basal"], 0.0)
        self.assertGreater(final_state["transduced_ciliated"], 0.0)
        self.assertGreater(final_state["total_transduced"], 30.0)

    def test_regulatory_and_structure_anchors(self) -> None:
        """Verify clinical trial identifiers, structural PDB tags, and mandatory NLM attribution."""
        self.assertEqual(CLINICAL_TRIAL_4D710, "NCT05248230")
        self.assertEqual(CLINICAL_TRIAL_VX522, "NCT05668741")
        self.assertEqual(PDB_DEPHOSPHO_CFTR, "6MSM")
        self.assertEqual(PDB_ACTIVE_CFTR, "5UAK")
        self.assertIn("ClinicalTrials.gov does not evaluate safety or efficacy", NLM_ATTRIBUTION)


if __name__ == "__main__":
    unittest.main()
