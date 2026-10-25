#!/usr/bin/env python3
"""Unit test suite for Airway Basal Stem Cell Quiescence Exit & Retrograde Lineage Kinetics (T-0414)."""

import math
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from basal_stem_quiescence_kinetics import (
    CLINICAL_TRIAL_4D710,
    CURATIVE_CONDUCTANCE_THRESHOLD_PCT,
    DEFAULT_EPISOMAL_HALFLIFE_DAYS,
    NLM_ATTRIBUTION,
    PDB_ACTIVE_CFTR,
    PRIMARY_DOI_BIOENG,
    PRIMARY_DOI_KIM,
    PRIMARY_DOI_SUK,
    BasalStemQuiescenceModel,
    run_full_simulation,
)


class TestBasalStemQuiescenceKinetics(unittest.TestCase):
    """Test suite verifying biophysical dynamics, lineage ODEs, and modality longevity."""

    def setUp(self) -> None:
        self.model = BasalStemQuiescenceModel(
            b0=30.0,
            quiescent_fraction=0.85,
            lambda_div=0.010,
            p_sym=0.25,
            episomal_halflife_days=100.0,
            curative_threshold=10.0,
        )

    def test_homeostatic_initialization_and_steady_state(self) -> None:
        """Verify homeostatic compartment sizing and near-zero baseline derivatives."""
        self.assertAlmostEqual(self.model.b0, 30.0, places=4)
        self.assertAlmostEqual(self.model.bq0, 25.5, places=4)
        self.assertAlmostEqual(self.model.ba0, 4.5, places=4)
        self.assertGreater(self.model.ta0, 0.0)
        self.assertGreater(self.model.diff0, 0.0)

        # Baseline derivatives at zero injury
        dbq, dba, dta, ddiff = self.model.derivatives(
            self.model.bq0, self.model.ba0, self.model.ta0, self.model.diff0, injury_level=0.0
        )
        self.assertAlmostEqual(dbq, 0.0, places=4)
        self.assertAlmostEqual(dta, 0.0, places=4)
        self.assertAlmostEqual(ddiff, 0.0, places=4)

        # Run 365-day uninjured simulation
        traj = self.model.simulate_hierarchy_rk4(days=365.0, dt=1.0, injury_level=0.0)
        self.assertEqual(len(traj), 366)
        final_t, final_bq, final_ba, final_ta, final_diff = traj[-1]
        final_b_total = final_bq + final_ba
        self.assertAlmostEqual(final_b_total, 30.0, delta=1.5)

    def test_injury_induced_quiescence_exit(self) -> None:
        """Verify injury drives quiescent stem cells into active cycling and expands progenitors."""
        dbq_base, dba_base, _, _ = self.model.derivatives(
            self.model.bq0, self.model.ba0, self.model.ta0, self.model.diff0, injury_level=0.0
        )
        dbq_inj, dba_inj, dta_inj, _ = self.model.derivatives(
            self.model.bq0, self.model.ba0, self.model.ta0, self.model.diff0, injury_level=1.0
        )

        # Under injury, quiescent stem cells exit (dbq becomes strongly negative, dba increases)
        self.assertLess(dbq_inj, dbq_base)
        self.assertGreater(dba_inj, dba_base)
        self.assertGreater(dta_inj, 0.0)

        # Dynamic simulation under transient injury shows recovery back toward homeostasis
        traj = self.model.simulate_hierarchy_rk4(days=180.0, dt=0.5, injury_level=2.0, injury_decay_days=10.0)
        t_final, bq_final, ba_final, _, _ = traj[-1]
        b_total_final = bq_final + ba_final
        self.assertAlmostEqual(b_total_final, 30.0, delta=3.0)

    def test_episomal_vector_washout_kinetics(self) -> None:
        """Verify exponential washout of non-integrating vectors matches t_1/2 ~ 100 days."""
        c0 = 50.0  # 50% initial transduction
        _, lum_0, cond_0 = self.model.calculate_transgene_retention(0.0, False, 10.0, c0)
        _, lum_100, cond_100 = self.model.calculate_transgene_retention(100.0, False, 10.0, c0)
        _, lum_200, cond_200 = self.model.calculate_transgene_retention(200.0, False, 10.0, c0)
        _, lum_1000, cond_1000 = self.model.calculate_transgene_retention(1000.0, False, 10.0, c0)

        self.assertAlmostEqual(lum_0, c0, places=4)
        self.assertAlmostEqual(cond_0, c0, places=4)
        self.assertAlmostEqual(lum_100, 0.5 * c0, places=3)
        self.assertAlmostEqual(cond_100, 0.5 * c0, places=3)
        self.assertAlmostEqual(cond_200, 0.25 * c0, places=3)
        self.assertLess(cond_1000, 0.1)

    def test_integrative_basal_editing_persistence(self) -> None:
        """Verify permanent basal integration reaches homeostatic asymptotic floor across 10 years."""
        basal_edited = 25.0
        luminal_init = 10.0

        b_0, lum_0, cond_0 = self.model.calculate_transgene_retention(0.0, True, basal_edited, luminal_init)
        b_365, lum_365, cond_365 = self.model.calculate_transgene_retention(365.0, True, basal_edited, luminal_init)
        b_3650, lum_3650, cond_3650 = self.model.calculate_transgene_retention(3650.0, True, basal_edited, luminal_init)

        # Basal stem cell editing remains exactly constant over decades
        self.assertAlmostEqual(b_0, basal_edited, places=4)
        self.assertAlmostEqual(b_365, basal_edited, places=4)
        self.assertAlmostEqual(b_3650, basal_edited, places=4)

        # Luminal epithelium starts at 10% and approaches 25% asymptotically
        self.assertAlmostEqual(lum_0, 10.0, places=4)
        self.assertGreater(lum_365, 23.0)
        self.assertAlmostEqual(lum_3650, 25.0, places=4)
        self.assertAlmostEqual(cond_3650, 25.0, places=4)
        self.assertGreaterEqual(cond_3650, CURATIVE_CONDUCTANCE_THRESHOLD_PCT)

    def test_redosing_interval_calculation(self) -> None:
        """Verify redosing interval derivation for non-integrating therapies."""
        # 40% initial, threshold 10%: two half-lives (40 -> 20 -> 10) = 200 days at t_1/2 = 100 days
        interval = self.model.calculate_redosing_interval_days(40.0, threshold_pct=10.0)
        self.assertIsNotNone(interval)
        self.assertAlmostEqual(interval, 200.0, delta=0.5)

        # Initial below threshold returns 0.0
        self.assertEqual(self.model.calculate_redosing_interval_days(8.0, threshold_pct=10.0), 0.0)

    def test_modality_comparisons_and_registry_grounding(self) -> None:
        """Verify modality benchmarking results and primary registry/literature grounding."""
        comparisons = self.model.compare_modalities(days=3650.0)
        self.assertEqual(len(comparisons), 3)

        mrna_mod, aav_mod, int_mod = comparisons

        # Episomal modalities wash out
        self.assertFalse(mrna_mod.is_integrative)
        self.assertFalse(aav_mod.is_integrative)
        self.assertLess(mrna_mod.year_10_conductance_pct_wt, 0.001)
        self.assertLess(aav_mod.year_10_conductance_pct_wt, 0.001)
        self.assertIsNotNone(mrna_mod.redosing_interval_days)
        self.assertIsNotNone(aav_mod.redosing_interval_days)

        # Integrative modality persists
        self.assertTrue(int_mod.is_integrative)
        self.assertAlmostEqual(int_mod.year_10_conductance_pct_wt, 25.0, places=1)
        self.assertIsNone(int_mod.redosing_interval_days)

    def test_full_simulation_schema_and_anti_hallucination_constants(self) -> None:
        """Verify full simulation dictionary contains all mandatory audit fields and citations."""
        res = run_full_simulation(days=3650.0)
        self.assertEqual(res["task_id"], "T-0414")
        self.assertEqual(res["epistemic_tier"], "Tier T2")
        self.assertIn("NLM", res["nlm_attribution"])
        self.assertIn("ClinicalTrials.gov", res["nlm_attribution"])
        self.assertIn(PRIMARY_DOI_SUK, res["primary_grounding"].values())
        self.assertIn(PRIMARY_DOI_BIOENG, res["primary_grounding"].values())
        self.assertIn(PRIMARY_DOI_KIM, res["primary_grounding"].values())
        self.assertEqual(res["primary_grounding"]["clinical_trial_nct"], CLINICAL_TRIAL_4D710)
        self.assertEqual(res["primary_grounding"]["pdb_structure"], PDB_ACTIVE_CFTR)

    def test_parameter_validation_bounds(self) -> None:
        """Verify constructor bounds checking for biological realism."""
        with self.assertRaises(ValueError):
            BasalStemQuiescenceModel(b0=-5.0)
        with self.assertRaises(ValueError):
            BasalStemQuiescenceModel(quiescent_fraction=1.5)
        with self.assertRaises(ValueError):
            BasalStemQuiescenceModel(lambda_div=0.50)  # unrealistic hyperproliferation
        with self.assertRaises(ValueError):
            BasalStemQuiescenceModel(p_sym=0.01)       # unviable asymmetric loss
        with self.assertRaises(ValueError):
            BasalStemQuiescenceModel(episomal_halflife_days=-10.0)


if __name__ == "__main__":
    unittest.main()
